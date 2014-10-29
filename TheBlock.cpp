#include "FinalSuperblock.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixX_t& hS,
                   const std::vector<MatrixX_t>& rhoBasisH2)
                   : m(m), hS(hS), rhoBasisH2(rhoBasisH2) {};

TheBlock::TheBlock(const Hamiltonian& ham) : m(d), hS(ham.h1)
{
    rhoBasisH2.assign(ham.siteBasisH2.begin(),
                      ham.siteBasisH2.begin() + nIndepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                             double& cumulativeTruncationError)
{
    MatrixX_t hSprime = createHprime(this, data.ham);  // expanded system block
    int md = m * d;
    if(data.exactDiag)
        return TheBlock(md, hSprime, createNewRhoBasisH2(data.ham.siteBasisH2,
                                                         true));
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    solveHSuper(hSprime, data, psiGround);            // calculate ground state
    int compm = data.compBlock -> m;
    psiGround.resize(md, compm * d);
    SelfAdjointEigenSolver<MatrixX_t> rhoSolver(psiGround * psiGround.adjoint());
                                             // find density matrix eigenstates
    int evecsToKeep;
    if(md <= data.mMax)
        evecsToKeep = md;
    else
    {
        int firstKeptEval = md - data.mMax;
        for(double currentFirstEval = rhoSolver.eigenvalues()(firstKeptEval);
            firstKeptEval < md
            && (currentFirstEval == 0
                || (currentFirstEval - rhoSolver.eigenvalues()(firstKeptEval - 1))
                   / std::abs(currentFirstEval) < degenerateDMCutoff);
            firstKeptEval++);
              // find the the max number of eigenvectors to keep that do not
              // terminate inside a degenerate eigenspace of the density matrix
        evecsToKeep = md - firstKeptEval;
        if(evecsToKeep == 0)
        {
            std::cerr << "More than mMax highest-weighted density-matrix "
                      << "eigenvectors are degenerate." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if(evecsToKeep != data.mMax)
            std::cout << "Warning: mMax truncation ends in a degenerate DM "
                      << "eigenspace, lowering cutoff to " << evecsToKeep
                      << " states." << std::endl;
    };
    cumulativeTruncationError
        += rhoSolver.eigenvalues().head(md - evecsToKeep).sum();
    primeToRhoBasis = rhoSolver.eigenvectors().rightCols(evecsToKeep);
                                            // construct change-of-basis matrix
    if(!data.infiniteStage) // modify psiGround to predict the next ground state
    {
        for(int sPrimeIndex = 0; sPrimeIndex < md; sPrimeIndex++)
                // transpose the environment block and the right-hand free site
        {
            rmMatrixX_t ePrime = psiGround.row(sPrimeIndex);
            ePrime.resize(compm, d);
            ePrime.transposeInPlace();
            ePrime.resize(1, d * compm);
            psiGround.row(sPrimeIndex) = ePrime;
        };
        psiGround = primeToRhoBasis.adjoint() * psiGround; 
                                      // change the expanded system block basis
        psiGround.resize(evecsToKeep * d, compm);
        psiGround *= data.beforeCompBlock -> primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(evecsToKeep * d * data.beforeCompBlock -> m * d, 1);
    };
    return TheBlock(evecsToKeep, changeBasis(hSprime),
                    createNewRhoBasisH2(data.ham.siteBasisH2, false));
                                  // save expanded-block operators in new basis
};

MatrixX_t TheBlock::createHprime(const TheBlock* block, const Hamiltonian& ham)
    const
{
    return kp(block -> hS, Id_d)
           + ham.blockSiteJoin(block -> rhoBasisH2)
           + kp(Id(block -> m), ham.h1);
};

std::vector<MatrixX_t> TheBlock::createNewRhoBasisH2(const vecMatD_t& siteBasisH2,
                                                     bool exactDiag) const
{
    std::vector<MatrixX_t> newRhoBasisH2;
    newRhoBasisH2.reserve(nIndepCouplingOperators);
    for(auto op = siteBasisH2.begin(), end = op + nIndepCouplingOperators;
        op != end; op++)
        newRhoBasisH2.push_back(exactDiag ?
                                kp(Id(m), *op) :
                                changeBasis(kp(Id(m), *op)));
    return newRhoBasisH2;
};

double TheBlock::solveHSuper(const MatrixX_t& hSprime, const stepData& data,
                             rmMatrixX_t& psiGround) const
{
    MatrixX_t hEprime = (data.infiniteStage ?
                         hSprime :
                         createHprime(data.compBlock, data.ham));
                                                  // expanded environment block
    MatrixX_t hSuper = kp(hSprime, Id(data.compBlock -> m * d))
                       + data.ham.siteSiteJoin(m, data.compBlock -> m)
                       + kp(Id(m * d), hEprime);
    return lanczos(hSuper, psiGround, data.lancTolerance);
};

MatrixX_t TheBlock::changeBasis(const MatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            rmMatrixX_t& psiGround, int lSys,
                                            int skips) const
{
    MatrixX_t hSprime = createHprime(this, data.ham);  // expanded system block
    double gsEnergy = solveHSuper(hSprime, data, psiGround);
                                                      // calculate ground state
    return FinalSuperblock(gsEnergy, lSys, psiGround, m, data.compBlock -> m,
                           skips);
};

obsMatrixX_t TheBlock::obsChangeBasis(const obsMatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
