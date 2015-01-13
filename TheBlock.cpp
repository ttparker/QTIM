#include "FinalSuperblock.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

TheBlock::TheBlock(const Hamiltonian& ham)
{
    blockParts.m = d;
    blockParts.hS = ham.h1;
    blockParts.rhoBasisH2
        .assign(ham.siteBasisH2.begin(),
                ham.siteBasisH2.begin() + nIndepCouplingOperators);
};

TheBlock::TheBlock(const effectiveHams& blockParts) : blockParts(blockParts) {};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                             double& cumulativeTruncationError)
{
    int md = blockParts.m * d;
    if(data.exactDiag)
    {
        effectiveHams newBlockParts;
        newBlockParts.m = md;
        newBlockParts.hS =   kp(blockParts.hS, Id_d)
                           + data.ham.EDNNCoupling(blockParts.rhoBasisH2)
                           + kp(MatrixXd::Identity(blockParts.m, blockParts.m),
                                data.ham.h1);
        newBlockParts.rhoBasisH2.reserve(nIndepCouplingOperators);
        for(auto op = data.ham.siteBasisH2.begin(), end = op + nIndepCouplingOperators;
            op != end; op++)
            newBlockParts.rhoBasisH2
                .push_back(kp(MatrixXd::Identity(blockParts.m, blockParts.m),
                              *op));
        return TheBlock(newBlockParts);
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    };
    lanczos(data.ham,
            data.infiniteStage ? blockParts : data.compBlock -> blockParts,
            psiGround, data.lancTolerance);           // calculate ground state
    int compm = data.compBlock -> blockParts.m;
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
        psiGround.resize(evecsToKeep * d
                         * data.beforeCompBlock -> blockParts.m * d, 1);
    };
    effectiveHams newBlockParts;
    newBlockParts.m = evecsToKeep;
    
    // project expanded system block down to highest-weighted DM eigenbasis:
    primeToRhoBasis.resize(blockParts.m, d * evecsToKeep);
    rmMatrixX_t oSO = blockParts.hS * primeToRhoBasis;
    oSO.resize(blockParts.m * d, evecsToKeep);
    primeToRhoBasis.resize(blockParts.m * d, evecsToKeep);
    newBlockParts.hS.noalias() = primeToRhoBasis.adjoint() * oSO;
                                                           // system block term
    newBlockParts.hS.noalias() += data.ham.NNCoupling(this);   // coupling term
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(evecsToKeep * blockParts.m, d);
    rmMatrixX_t oDagH1 = oDag * data.ham.h1;
    oDagH1.resize(evecsToKeep, blockParts.m * d);
    newBlockParts.hS.noalias() += oDagH1 * primeToRhoBasis;   // free site term
    
    // project coupling operators into new DM eigenspace basis:
    newBlockParts.rhoBasisH2.reserve(nIndepCouplingOperators);
    for(auto op = data.ham.siteBasisH2.begin(),
             end = op + nIndepCouplingOperators; op != end; op++)
    {
        rmMatrixX_t reshapedODag = primeToRhoBasis.adjoint();
        reshapedODag.resize(evecsToKeep * blockParts.m, d);
        rmMatrixX_t oDagH2 = reshapedODag * *op;
        oDagH2.resize(evecsToKeep, blockParts.m * d);
        newBlockParts.rhoBasisH2.push_back(oDagH2 * primeToRhoBasis);
    };
    return TheBlock(newBlockParts); // save expanded-block operators in new basis
};

rmMatrixX_t TheBlock::projectNNCoupling(const rmMatrixX_t& blockOp,
                                        const rmMatrixX_t& siteOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    rmMatrixX_t oDagHSite = oDag * siteOp,
                hSysO = blockOp * primeToRhoBasis;
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    oDagHSite.resize(nextSiteM, blockParts.m * d);
    hSysO.resize(blockParts.m * d, nextSiteM);
    return oDagHSite * hSysO;
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            rmMatrixX_t& psiGround, int lSys,
                                            int skips) const
{
    double gsEnergy = lanczos(data.ham, data.compBlock -> blockParts, psiGround,
                              data.lancTolerance);    // calculate ground state
    return FinalSuperblock(gsEnergy, lSys, psiGround, blockParts.m,
                           data.compBlock -> blockParts.m, skips);
};

obsMatrixX_t TheBlock::obsProjectBlockOp(const obsMatrixX_t& sysOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    obsMatrixX_t oSO = sysOp * primeToRhoBasis;
    oSO.resize(blockParts.m * d, nextSiteM);
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    return primeToRhoBasis.adjoint() * oSO;
};

#ifdef differentScalars
obsMatrixX_t TheBlock::obsProjectBlockFreeSiteOps(const obsMatrixX_t& blockOp,
                                                  const obsMatrixX_t& siteOp)
{
    int nextSiteM = primeToRhoBasis.cols();
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    primeToRhoBasis.resize(blockParts.m, d * nextSiteM);
    obsMatrixX_t oDagHSite = oDag * siteOp,
                 hSysO = blockOp * primeToRhoBasis;
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    oDagHSite.resize(nextSiteM, blockParts.m * d);
    hSysO.resize(blockParts.m * d, nextSiteM);
    return oDagHSite * hSysO;
};
#endif

obsMatrixX_t TheBlock::obsProjectFreeSiteOp(const obsMatrixD_t& lFreeSite)
{
    int nextSiteM = primeToRhoBasis.cols();
    primeToRhoBasis.resize(blockParts.m * d, nextSiteM);
    rmMatrixX_t oDag = primeToRhoBasis.adjoint();
    oDag.resize(nextSiteM * blockParts.m, d);
    obsMatrixX_t oDagH1 = oDag * lFreeSite;
    oDagH1.resize(nextSiteM, blockParts.m * d);
    return oDagH1 * primeToRhoBasis;
};
