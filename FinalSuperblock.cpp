#include "FinalSuperblock.h"
#include "Lanczos.h"

using namespace Eigen;

FinalSuperblock::FinalSuperblock(const MatrixXd& matFinal,
                                 const rmMatrixXd& psiGroundIn,
                                 stepData data, int mSFinal, int mEFinal,
                                 int skips)
    : lSupFinal(data.ham.lSys), psiGround(psiGroundIn), mSFinal(mSFinal),
      mEFinal(mEFinal), skips(skips)
{
    gsEnergy = lanczos(matFinal, psiGround, data.lancTolerance);
    if(lSupFinal % 2)
    {
        lSFinal = (lSupFinal - 1)/2;
        lEFinal = (lSupFinal - 3)/2;
    }
    else
        lSFinal = lEFinal = lSupFinal / 2 - 1;
};

double FinalSuperblock::expValue(const opsVec& ops,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks)
{
    opsMap sysBlockOps, // observable operators that will act on the system block
           envBlockOps; // same for environment block
    MatrixDd lFreeSite = Id_d,
             rFreeSite = Id_d;
    bool opAtlFreeSite = false,
         opAtrFreeSite = false;
    for(const auto& opToEvaluate : ops) // split up operators into their blocks
    {
        int site = opToEvaluate.second;
        if(site < lSFinal)
            placeOp(opToEvaluate, sysBlockOps, true);
        else if(site == lSFinal)
        {                       // if necessary, assign left free-site operator
            lFreeSite *= opToEvaluate.first;
            opAtlFreeSite = true;
        }
        else if(site == lSFinal + 1)
        {
            rFreeSite *= opToEvaluate.first;    // and right free-site operator
            opAtrFreeSite = true;
        }
        else
            placeOp(opToEvaluate, envBlockOps, false);
    };
    if(sysBlockOps.empty() && !opAtlFreeSite)
                        // observables in right-hand half of superblock case
        if(envBlockOps.empty())
        {                   // right-hand free site single-site observable case
            psiGround.resize(mSFinal * d * mEFinal, d);
            return (psiGround
                    * rFreeSite.transpose()
                    * psiGround.adjoint()
                   ).trace();
        }
        else
        {
            psiGround.resize(mSFinal * d, mEFinal * d);
            return (psiGround.adjoint()
                    * psiGround
                    * kp(rhoBasisRep(envBlockOps, rightBlocks, lEFinal),
                         rFreeSite).transpose()
                   ).trace();
        }
    else if(envBlockOps.empty() && !opAtrFreeSite)
                            // observables in left-hand half of superblock case
        if(!opAtlFreeSite)              // all observables in system block case
        {
            psiGround.resize(mSFinal, d * mEFinal * d);
            return (psiGround.adjoint()
                    * rhoBasisRep(sysBlockOps, leftBlocks, lSFinal)
                    * psiGround
                   ).trace();
        }
        else
        {
            psiGround.resize(mSFinal * d, mEFinal * d);
            return (psiGround.adjoint()
                    * kp(rhoBasisRep(sysBlockOps, leftBlocks, lSFinal),
                         lFreeSite)
                    * psiGround
                   ).trace();
        }
    else                    // observables in both halves of superblock case
    {
        psiGround.resize(mSFinal * d, mEFinal * d);
        return (psiGround.adjoint()
                * kp(rhoBasisRep(sysBlockOps, leftBlocks, lSFinal), lFreeSite)
                * psiGround
                * kp(rhoBasisRep(envBlockOps, rightBlocks, lEFinal), rFreeSite)
                  .transpose()
               ).trace();
    };
};

void FinalSuperblock::placeOp(const std::pair<MatrixDd, int>& op,
                              opsMap& blockSide, bool systemSide)
{
    int lhSite = (systemSide ? op.second : lSupFinal - 1 - op.second);
    if(blockSide.count(lhSite))      // already an observable at this site?
        blockSide[lhSite] *= op.first;
    else
        blockSide.insert(std::pair<int, MatrixDd>(lhSite, op.first));
};

MatrixXd FinalSuperblock::rhoBasisRep(const opsMap& blockOps,
                                      std::vector<TheBlock>& blocks,
                                      int blockSize) const
{
    if(blockOps.empty())
        return Id(blocks[blockSize - 1].m);
    MatrixXd rhoBasisBlockOp;
    const auto firstOp = blockOps.begin();
                // first nontrivial site operator at which to start tensoring
    auto currentOp = firstOp;
    const auto firstSite = blocks.begin();
    auto currentSite = firstSite;
    if(firstOp -> first == 0)       // one-site observable at first site?
    {
        rhoBasisBlockOp = firstOp -> second;
        currentOp++;
    }
    else
    {
        currentSite += (firstOp -> first - 1);
        rhoBasisBlockOp = Id(currentSite -> m);
    };
    const auto opsEnd = blockOps.end();
    for(const auto end = firstSite + (blockSize - 1); currentSite != end;
        currentSite++)
    {
        MatrixXd primeBasisBlockOp = kp(rhoBasisBlockOp,
                                        (currentOp != opsEnd
                                         && currentSite - firstSite
                                            == currentOp -> first - 1 ?
                                         currentOp++ -> second :
                                         Id_d));
        rhoBasisBlockOp = (currentSite - firstSite < skips ?
                           primeBasisBlockOp :
                           currentSite -> changeBasis(primeBasisBlockOp));
    };
    return rhoBasisBlockOp;
};