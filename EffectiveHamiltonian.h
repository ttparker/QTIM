#ifndef EFFECTIVEHAMILTONAIN_H
#define EFFECTIVEHAMILTONAIN_H

#include <map>
#include "TheBlock.h"

typedef std::vector<std::pair<MatrixDd, int>,
                    Eigen::aligned_allocator<std::pair<MatrixDd, int>>> opsVec;
typedef std::map<int, MatrixDd, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, MatrixDd>>> opsMap;

class EffectiveHamiltonian
{
    public:
        EffectiveHamiltonian(const Eigen::MatrixXd& matFinal, int lSupFinal,
                             int mSFinal, int skips);
        double gsEnergy() const,                 // returns ground-state energy
               expValue(const opsVec& ops, std::vector<TheBlock>& leftBlocks,
                        std::vector<TheBlock>& rightBlocks);
        // calculates exectation value of a combination of single-site operators

    private:
        double storedGSEnergy;
        rmMatrixXd psiGround;				// final superblock ground state
        int lSupFinal,                      // final system size
            mSFinal,                // final number of states stored per block
            skips;
        
        void placeOp(const std::pair<MatrixDd, int>& op, opsMap& blockSide,
                     bool systemSide);
                    // assign each one-site observable to the appropriate block
        Eigen::MatrixXd rhoBasisRep(const opsMap& blockOps,
                                    std::vector<TheBlock>& blocks) const;
                // converts single-site operators into the system block basis
};

#endif
