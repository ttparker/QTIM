#include "FinalSuperblock.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixX_t& hS,
                   const std::vector<MatrixX_t>& rhoBasisH2)
                   : m(m), hS(hS), rhoBasisH2(rhoBasisH2) {};

TheBlock::TheBlock(const Hamiltonian& ham, bool westSide)
    : m(d), hS(westSide ? ham.westSideH1 : ham.eastSideH1)
{
    rhoBasisH2.assign(ham.siteBasisH2.begin(),
                      ham.siteBasisH2.begin() + indepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                             TheBlock* nextCompBlock)
{
    MatrixX_t hSprime = createHprime(this, true, data), // expanded system block
              hEprime = createHprime(data.compBlock, false, data);
                                                  // expanded environment block
    int md = m * d;
    if(data.exactDiag)
    {
        *nextCompBlock = TheBlock(md, hEprime,
                                  createNewRhoBasisH2(data.ham.siteBasisH2,
                                                      true, this));
        return TheBlock(md, hSprime, createNewRhoBasisH2(data.ham.siteBasisH2,
                                                         true, this));
    };
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    solveHSuper(hSprime, hEprime, data, psiGround);   // calculate ground state
    int compm = data.compBlock -> m;
    psiGround.resize(md, compm * d);
    SelfAdjointEigenSolver<MatrixX_t> rhoSolver(psiGround * psiGround.adjoint());
                                      // find system density matrix eigenstates
    primeToRhoBasis = rhoSolver.eigenvectors().rightCols(data.mMax);
                                     // construct system change-of-basis matrix
    if(data.infiniteStage)
    {
        rhoSolver.compute(psiGround.adjoint() * psiGround);
                                 // find environment density matrix eigenstates
        data.compBlock -> primeToRhoBasis
            = rhoSolver.eigenvectors().rightCols(data.mMax);
                                // construct environment change-of-basis matrix
        *nextCompBlock = TheBlock(data.mMax,
                                  changeBasis(hEprime, data.compBlock),
                                  createNewRhoBasisH2(data.ham.siteBasisH2,
                                                      false, data.compBlock));
    }
    else                   // modify psiGround to predict the next ground state
    {
        for(int sPrimeIndex = 0; sPrimeIndex < md; sPrimeIndex++)
                    // transpose the environment block and right-hand free site
        {
            rmMatrixX_t ePrime = psiGround.row(sPrimeIndex);
            ePrime.resize(compm, d);
            ePrime.transposeInPlace();
            ePrime.resize(1, compm * d);
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
    return TheBlock(data.mMax, changeBasis(hSprime, this),
                    createNewRhoBasisH2(data.ham.siteBasisH2, false, this));
                                  // save expanded-block operators in new basis
};

MatrixX_t TheBlock::createHprime(const TheBlock* block, bool hSprime,
                                 const stepData& data) const
{
    return kp(block -> hS, Id_d)
           + data.ham.blockSiteJoin(block -> rhoBasisH2)
           + kp(Id(block -> m), hSprime == data.sweepingEast ?
                                data.ham.westSideH1 :
                                data.ham.eastSideH1);
};

std::vector<MatrixX_t> TheBlock::createNewRhoBasisH2(const vecMatD_t& siteBasisH2,
                                                     bool exactDiag,
                                                     const TheBlock* block) const
{
    std::vector<MatrixX_t> newRhoBasisH2;
    newRhoBasisH2.reserve(indepCouplingOperators);
    for(auto op = siteBasisH2.begin(), end = op + indepCouplingOperators;
        op != end; op++)
        newRhoBasisH2.push_back(exactDiag ?
                                kp(Id(m), *op) :
                                changeBasis(kp(Id(m), *op), block));
    return newRhoBasisH2;
};

double TheBlock::solveHSuper(const MatrixX_t& hSprime, const MatrixX_t& hEprime,
                             const stepData& data, rmMatrixX_t& psiGround) const
{
    MatrixX_t hSuper = kp(hSprime, Id(data.compBlock -> m * d))
                       + data.ham.siteSiteJoin(m, data.compBlock -> m)
                       + kp(Id(m * d), hEprime);
    return lanczos(hSuper, psiGround, data.lancTolerance);
};

MatrixX_t TheBlock::changeBasis(const MatrixX_t& mat, const TheBlock* block)
                                const
{
    return block -> primeToRhoBasis.adjoint() * mat * block -> primeToRhoBasis;
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            rmMatrixX_t& psiGround, int skips)
                                            const
{
    MatrixX_t hSprime = createHprime(this, true, data), // expanded system block
              hEprime = createHprime(data.compBlock, false, data);
                                                  // expanded environment block
    double gsEnergy = solveHSuper(hSprime, hEprime, data, psiGround);
                                                      // calculate ground state
    return FinalSuperblock(gsEnergy, data.ham.lSys, psiGround, m,
                           data.compBlock -> m, skips);
};

obsMatrixX_t TheBlock::obsChangeBasis(const obsMatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
