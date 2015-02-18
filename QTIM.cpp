#include "TheBlock.h"

#define j couplingConstants[0]
#define h couplingConstantsIn[1]
#define rhoBasisSigmax rhoBasisH2[0]
#define siteBasisSigmax siteBasisH2[0]

#define mS blockParts.m
#define mE compBlockParts.m
#define hSiteBond siteBasisH2[0]
#define hSBond blockParts.rhoBasisH2[0]
#define hEBond compBlockParts.rhoBasisH2[0]

using namespace Eigen;

Hamiltonian::Hamiltonian()
{
    siteBasisH2.resize(nCouplingOperators);
    hSiteBond << 0., .5,
                 .5, 0.;
};

void Hamiltonian::setParams(const std::vector<double>& couplingConstantsIn)
{
    h1 << -h / 2,    0.,
              0., h / 2;
    couplingConstants = {couplingConstantsIn[0]};
};

rmMatrixX_t Hamiltonian::EDNNCoupling(const std::vector<MatrixX_t>& rhoBasisH2)
    const
{
    return j * kp(rhoBasisSigmax, siteBasisSigmax);
};

Matrix<hamScalar, Dynamic, 1>
    Hamiltonian::act(const effectiveHams& blockParts,
                     const effectiveHams& compBlockParts, rmMatrixX_t x) const
{
    int mSd = mS * d,
        mEd = mE * d,
        totalDimension = mSd * mEd;
                                // dimension of state on which Hamiltonian acts
    
    // act with right-hand half of system:
    x.resize(mSd * mE, d);
    rmMatrixX_t cycledx = x.transpose(),
                y4123 = h1 * cycledx,          // act with right-hand free site
                hBondx = hSiteBond * cycledx;
                                 // act with right-hand free site half of
                                 // environment-block-right-hand-free-site bond
    cycledx.resize(d * mSd, mE);
    hBondx.resize(d * mSd, mE);
    y4123.resize(d * mSd, mE);
    y4123.noalias() +=   cycledx * compBlockParts.hS.transpose()
                                              // act with environment block ...
                       + j * hBondx * hEBond.transpose();
                                 // ... and with environment-block half of
                                 // environment-block-right-hand-free-site bond
    y4123.resize(d, mSd * mE);
    y4123.transposeInPlace();
    y4123.resize(totalDimension, 1);
    
    // act with left-hand side of system:
    x.resize(mS, d * mEd);
    rmMatrixX_t y2341 = (blockParts.hS * x).transpose(); // act with system block
    cycledx = x.transpose();
    hBondx = cycledx * hSBond.transpose();
         // act with system-block half of system-block-left-hand-free-site bond
    cycledx.resize(d, mEd * mS);
    hBondx.resize(d, mEd * mS);
    y2341.resize(d, mEd * mS);
    y2341.noalias() +=   h1 * cycledx   // act with left-hand free site ...
                       + j * hSiteBond * hBondx;
                                    // ... and with left-hand-free-site half of
                                    // system-block-left-hand-free-site bond
    y2341.resize(d * mEd, mS);
    y2341.transposeInPlace();
    y2341.resize(totalDimension, 1);
    
    // act with free sites' bond:
    rmMatrixX_t yFSBond(mS, d * mEd);
    for(int i = 0; i < mS; i++)
    {
        rmMatrixX_t RHS = x.row(i);
        RHS.resize(d, mEd);
        rmMatrixX_t hSiteBondRHS = hSiteBond * RHS;
                // act with left-hand-free-site half of bond between free sites
        hSiteBondRHS.resize(mEd, d);
        rmMatrixX_t hSiteBondRHShSiteBond
            = j * hSiteBondRHS * hSiteBond.transpose();
                // act with left-hand-free-site half of bond between free sites
        hSiteBondRHShSiteBond.resize(1, d * mEd);
        yFSBond.row(i) = hSiteBondRHShSiteBond;
    };
    yFSBond.resize(totalDimension, 1);
    
    return y4123 + y2341 + yFSBond;
};

rmMatrixX_t Hamiltonian::NNCoupling(TheBlock * block) const
{
    return j * block -> projectNNCoupling(block -> blockParts.rhoBasisSigmax,
                                          siteBasisSigmax);
};
