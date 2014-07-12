#include "FinalSuperblock.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixX_t& hS,
                   const std::vector<MatrixX_t>& rhoBasisH2)
                   : m(m), hS(hS), rhoBasisH2(rhoBasisH2) {};

TheBlock::TheBlock(const Hamiltonian& ham) : m(d), hS(ham.h1)
{
    rhoBasisH2.assign(ham.siteBasisH2.begin(),
                      ham.siteBasisH2.begin() + indepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround)
{
    MatrixX_t hSprime = createHprime(this, data.ham);  // expanded system block
    int md = m * d;
    if(data.exactDiag)
        return TheBlock(md, hSprime, createNewRhoBasisH2(data.ham.siteBasisH2,
                                                         true));
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    int compm = data.compBlock -> m,
        compmd = compm * d;
    MatrixX_t hEprime = (data.infiniteStage ?
                         hSprime :
                         createHprime(data.compBlock, data.ham));
                                                  // expanded environment block
    lanczos(kp(hSprime, Id(compmd))
               + data.ham.siteSiteJoin(m, compm)
               + kp(Id(md), hEprime),
            psiGround, data.lancTolerance);	          // calculate ground state
    psiGround.resize(md, compmd);
    SelfAdjointEigenSolver<MatrixX_t> rhoSolver(psiGround * psiGround.adjoint());
                                             // find density matrix eigenstates
    primeToRhoBasis = rhoSolver.eigenvectors().rightCols(data.mMax);
                                            // construct change-of-basis matrix
    if(!data.infiniteStage) // modify psiGround to predict the next ground state
    {
        for(int sPrimeIndex = 0; sPrimeIndex < md; sPrimeIndex++)
                    // transpose the environment block and right-hand free site
        {
            rmMatrixX_t ePrime = psiGround.row(sPrimeIndex);
            ePrime.resize(compm, d);
            ePrime.transposeInPlace();
            ePrime.resize(1, compmd);
            psiGround.row(sPrimeIndex) = ePrime;
        };
        psiGround = primeToRhoBasis.adjoint() * psiGround; 
                                      // change the expanded system block basis
        psiGround.resize(data.mMax * d, compm);
        psiGround *= data.beforeCompBlock -> primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(data.mMax * d
                         * data.beforeCompBlock -> primeToRhoBasis.rows(), 1);
    };
    return TheBlock(data.mMax, changeBasis(hSprime),
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
    newRhoBasisH2.reserve(indepCouplingOperators);
    for(auto op = siteBasisH2.begin(), end = op + indepCouplingOperators;
        op != end; op++)
        newRhoBasisH2.push_back(exactDiag ?
                                kp(Id(m), *op) :
                                changeBasis(kp(Id(m), *op)));
    return newRhoBasisH2;
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            rmMatrixX_t& psiGround, int skips)
                                            const
{
    MatrixX_t hSprime = createHprime(this, data.ham);  // expanded system block
    int compm = data.compBlock -> m;
    MatrixX_t hEprime = createHprime(data.compBlock, data.ham);
                                                  // expanded environment block
    MatrixX_t hSuper = kp(hSprime, Id(compm * d))
                          + data.ham.siteSiteJoin(m, compm)
                          + kp(Id(m * d), hEprime);
    double gsEnergy = lanczos(hSuper, psiGround, data.lancTolerance);
    return FinalSuperblock(gsEnergy, data.ham.lSys, psiGround, m, compm, skips);
};

obsMatrixX_t TheBlock::obsChangeBasis(const obsMatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};

MatrixX_t TheBlock::changeBasis(const MatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
