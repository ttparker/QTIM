#ifndef FINALSUPERBLOCK_H
#define FINALSUPERBLOCK_H

#include <map>
#include "TheBlock.h"

typedef std::vector<std::pair<obsMatrixD_t, int>,
                    Eigen::aligned_allocator<std::pair<obsMatrixD_t, int>>>
                    opsVec;
typedef std::map<int, obsMatrixD_t, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, obsMatrixD_t>>>
                 opsMap;

class FinalSuperblock
{
    public:
        double gsEnergy;                         // returns ground-state energy
        
        FinalSuperblock(const Eigen::MatrixXd& matFinal,
                        const rmMatrixXd& psiGroundIn, stepData data,
                        int mSFinal, int mEFinal, int skips);
        double expValue(const opsVec& ops, std::vector<TheBlock>& leftBlocks,
                        std::vector<TheBlock>& rightBlocks);
        // calculates exectation value of a combination of single-site operators
    
    private:
        int lSupFinal;                                     // final system size
        rmMatrixXd psiGround;                  // final superblock ground state
        int lSFinal,                        // final length of the system block
            lEFinal,                   // final length of the environment block
            mSFinal,           // final number of states stored in system block
            mEFinal,      // final number of states stored in environment block
            skips;
        
        void placeOp(const std::pair<obsMatrixD_t, int>& op, opsMap& blockSide,
                     bool systemSide);
                    // assign each one-site observable to the appropriate block
        obsMatrixX_t rhoBasisRep(const opsMap& blockOps,
                                 std::vector<TheBlock>& blocks, int blockSize)
                                 const;
                // converts single-site operators into the system block basis
};

#endif
