#ifndef THEBLOCK_H
#define THEBLOCK_H

#include "Hamiltonian.h"
#define Id_d Matrix<double, d, d>::Identity()       // one-site identity matrix

class FinalSuperblock;
class TheBlock;

struct stepData
{
    Hamiltonian ham;                             // model Hamiltonian paramters
    bool exactDiag;             // close enough to edge to skip DMRG trucation?
    TheBlock* compBlock;         // complementary block on other side of system
    bool infiniteStage;
    double lancTolerance;  // max deviation from 1 of dot product of successive
                           // Lanczos iterations' ground state vectors
    int mMax;                                             // max bond dimension
    TheBlock* beforeCompBlock;     // next smaller block than complementary one
};

class TheBlock
{
    public:
        int m;                              // number of states stored in block
        MatrixX_t primeToRhoBasis;                    // change-of-basis matrix
        
        TheBlock(int m = 0, const MatrixX_t& hS = MatrixX_t(),
                 const std::vector<MatrixX_t>& rhoBasisH2
                     = std::vector<MatrixX_t>());
        TheBlock(const Hamiltonian& ham);
        TheBlock nextBlock(const stepData& data, rmMatrixX_t& psiGround);
                                                     // performs each DMRG step
        FinalSuperblock createHSuperFinal(const stepData& data,
                                          rmMatrixX_t& psiGround,
                                          int skips);
        obsMatrixX_t obsChangeBasis(const obsMatrixX_t& mat) const;
                       // changes basis during calculation of observables stage
    
    private:
        MatrixX_t hS;                                      // block Hamiltonian
        std::vector<MatrixX_t> rhoBasisH2;
                                     // density-matrix-basis coupling operators
        
        double lanczos(const MatrixX_t& mat, rmMatrixX_t& seed,
                       double lancTolerance);
     // changes input seed to ground eigenvector - make sure seed is normalized
        MatrixX_t createHprime(const TheBlock* block, const Hamiltonian& ham);
        MatrixX_t changeBasis(const MatrixX_t& mat) const;
                   // represents operators in the basis of the new system block
};

#endif
