#ifndef THEBLOCK_H
#define THEBLOCK_H

#include "Hamiltonian.h"
#define Id_d Matrix<double, d, d>::Identity()       // one-site identity matrix

class FinalSuperblock;

struct stepData
{
    Hamiltonian ham;                             // model Hamiltonian paramters
    bool exactDiag;             // close enough to edge to skip DMRG trucation?
    TheBlock* compBlock;         // complementary block on other side of system
    bool infiniteStage;
    double lancTolerance;  // max deviation from 1 of dot product of successive
                           // Lanczos iterations' ground state vectors
    int mMax;                              // max size of effective Hamiltonian
    TheBlock* beforeCompBlock;     // next smaller block than complementary one
};

class TheBlock
{
    public:
        effectiveHams blockParts; // stores effective Hamiltonians for this site
        rmMatrixX_t primeToRhoBasis;                  // change-of-basis matrix
        
        TheBlock() {};
        TheBlock(const Hamiltonian& ham);
        TheBlock(const effectiveHams& blockParts);
        TheBlock nextBlock(const stepData& data, rmMatrixX_t& psiGround,
                           double& cumulativeTruncationError);
                                                     // performs each DMRG step
        rmMatrixX_t projectNNCoupling(const rmMatrixX_t& blockOp,
                                      const rmMatrixX_t& siteOp);
                  // projects term coupling system block to left-hand free site
        FinalSuperblock createHSuperFinal(const stepData& data,
                                          rmMatrixX_t& psiGround,
                                          int lSys, int skips) const;
        obsMatrixX_t obsProjectBlockOp(const obsMatrixX_t& sysOp),
                      // projects any block op except for (the one at the inner
                      // end IF there is an op at the adjacent free site)
        #ifdef differentScalars
                     obsProjectNNCoupling(const obsMatrixX_t& blockOp,
                                          const obsMatrixX_t& siteOp),
                   // identical to member function above except for return type
        #endif
                     obsProjectFreeSiteOp(const obsMatrixD_t& lFreeSite);
                      // projects op at free site that isn't next to a block op
    
    private:
        double lanczos(const Hamiltonian& ham,
                       const effectiveHams& compBlockParts,
                       rmMatrixX_t& seed, double lancTolerance) const;
     // changes input seed to ground eigenvector - make sure seed is normalized
};

#endif
