#include "FinalSuperblock.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixX_t& hS, const matPair newRhoBasisH2s,
                   int l)
    : m(m), hS(hS), off0RhoBasisH2(newRhoBasisH2s.first),
      off1RhoBasisH2(newRhoBasisH2s.second), l(l) {};

TheBlock::TheBlock(const Hamiltonian& ham, bool westSide)
    : m(d), hS(westSide ? ham.westSideH1 : ham.eastSideH1), l(0)
{
    off0RhoBasisH2.assign(ham.siteBasisH2.begin(),
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
                                  createNewRhoBasisH2s(data.ham.siteBasisH2,
                                                      true, this), l + 1);
        return TheBlock(md, hSprime, createNewRhoBasisH2s(data.ham.siteBasisH2,
                                                         true, this), l + 1);
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
                                  createNewRhoBasisH2s(data.ham.siteBasisH2,
                                                       false, data.compBlock),
                                  l + 1);
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
                    createNewRhoBasisH2s(data.ham.siteBasisH2, false, this),
                    l + 1);       // save expanded-block operators in new basis
};

MatrixX_t TheBlock::createHprime(const TheBlock* block, bool hSprime,
                                 const stepData& data) const
{
    MatrixX_t hPrime = kp(block -> hS, Id_d)
                       + data.ham.blockAdjacentSiteJoin(1, block
                                                           -> off0RhoBasisH2)
                       + kp(Id(block -> m), hSprime == data.sweepingEast ?
                                            data.ham.westSideH1 :
                                            data.ham.eastSideH1);
    if(block -> l != 0)
        hPrime += data.ham.blockAdjacentSiteJoin(2, block -> off1RhoBasisH2);
    return hPrime;
};

matPair TheBlock::createNewRhoBasisH2s(const vecMatD_t& siteBasisH2,
                                       bool exactDiag, const TheBlock* block)
                                       const
{
    std::vector<MatrixX_t> newOff0RhoBasisH2,
                           newOff1RhoBasisH2;
    newOff0RhoBasisH2.reserve(indepCouplingOperators);
    newOff1RhoBasisH2.reserve(indepCouplingOperators);
    for(int i = 0; i < indepCouplingOperators; i++)
    {
        newOff0RhoBasisH2.push_back(exactDiag ?
                                    kp(Id(m), siteBasisH2[i]) :
                                    changeBasis(kp(Id(m), siteBasisH2[i]),
                                                block));
        newOff1RhoBasisH2.push_back(exactDiag ?
                                    kp(off0RhoBasisH2[i], Id_d) :
                                    changeBasis(kp(off0RhoBasisH2[i], Id_d),
                                                block));
    };
    return std::make_pair(newOff0RhoBasisH2, newOff1RhoBasisH2);
};

double TheBlock::solveHSuper(const MatrixX_t& hSprime, const MatrixX_t& hEprime,
                             const stepData& data, rmMatrixX_t& psiGround) const
{
    int compm = data.compBlock -> m;
    MatrixX_t hSuper = kp(hSprime, Id(compm * d))
                       + data.ham.lBlockrSiteJoin(off0RhoBasisH2, compm)
                       + data.ham.siteSiteJoin(m, compm)
                       + data.ham.lSiterBlockJoin(m, data.compBlock
                                                     -> off0RhoBasisH2)
                       + kp(Id(m * d), hEprime);                  // superblock
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
