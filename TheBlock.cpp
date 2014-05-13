#include "FinalSuperblock.h"
#include "Lanczos.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixX_t& hS,
                   const std::vector<MatrixX_t>& rhoBasisH2)
                   : m(m), hS(hS), rhoBasisH2(rhoBasisH2) {};

TheBlock::TheBlock(const Hamiltonian& ham) : m(d), hS(ham.h1)
{
    rhoBasisH2.assign(ham.h2.begin(), ham.h2.begin() + indepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround)
{
    MatrixX_t hSprime = kp(hS, Id_d)
                        + data.ham.blockSiteJoin(rhoBasisH2)
                        + kp(Id(m), data.ham.h1);      // expanded system block
    std::vector<MatrixX_t> tempRhoBasisH2;
    tempRhoBasisH2.reserve(indepCouplingOperators);
    int md = m * d;
    if(data.exactDiag)
    { // if near edge of system, no truncation necessary so skip DMRG algorithm
        for(auto op = data.ham.h2.begin(), end = op + indepCouplingOperators;
            op != end; op++)
            tempRhoBasisH2.push_back(kp(Id(m), *op));
        return TheBlock(md, hSprime, tempRhoBasisH2);
    };
    int compm = data.compBlock -> m,
        compmd = compm * d;
    lanczos(data.infiniteStage ?
            MatrixX_t(kp(hSprime, Id(md))
                      + data.ham.siteSiteJoin(m, m)
                      + kp(Id(md), hSprime)) :
            MatrixX_t(kp(hSprime, Id(compmd))
                      + data.ham.siteSiteJoin(m, compm)
                      + kp(Id(md), kp(Id(compm), data.ham.h1)
                                   + data.ham.blockSiteJoin(data.compBlock
                                                            -> rhoBasisH2)
                                   + kp(data.compBlock -> hS, Id_d))),
            psiGround, data.lancTolerance);	          // calculate ground state
    psiGround.resize(md, compmd);
    SelfAdjointEigenSolver<MatrixX_t> rhoSolver(psiGround * psiGround.adjoint());
                                             // find density matrix eigenstates
    primeToRhoBasis = rhoSolver.eigenvectors().rightCols(data.mMax);
                                            // construct change-of-basis matrix
    for(auto op = data.ham.h2.begin(), end = op + indepCouplingOperators;
        op != end; op++)
        tempRhoBasisH2.push_back(changeBasis(kp(Id(m), *op)));
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
    return TheBlock(data.mMax, changeBasis(hSprime), tempRhoBasisH2);
                                  // save expanded-block operators in new basis
};

obsMatrixX_t TheBlock::obsChangeBasis(const obsMatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            const rmMatrixX_t& psiGround,
                                            int skips) const
{
    int compm = data.compBlock -> m;
    return FinalSuperblock(kp(kp(hS, Id_d)
                              + data.ham.blockSiteJoin(rhoBasisH2)
                              + kp(Id(m), data.ham.h1),
                              Id(compm * d))
                           + data.ham.siteSiteJoin(m, compm)
                           + kp(Id(m * d), kp(Id(compm), data.ham.h1)
                                           + data.ham.blockSiteJoin(data.compBlock -> rhoBasisH2)
                                           + kp(data.compBlock -> hS, Id_d)),
                           psiGround, data, m, compm, skips);
};

MatrixX_t TheBlock::changeBasis(const MatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
