#ifndef THEBLOCK_H
#define THEBLOCK_H

#include "Hamiltonian.h"
#define Id_d MatrixDd::Identity()   // one-site identity matrix

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
    int mMax;                              // max size of effective Hamiltonian
    TheBlock* beforeCompBlock;     // next smaller block than complementary one
};

class TheBlock
{
    public:
        int m;                              // number of states stored in block
        Eigen::MatrixXd primeToRhoBasis;              // change-of-basis matrix
        
        TheBlock(int m = 0,
                 const Eigen::MatrixXd& hS = Eigen::MatrixXd(),
                 const std::vector<Eigen::MatrixXd>& rhoBasisH2
                     = std::vector<Eigen::MatrixXd>());
        TheBlock(const Hamiltonian& ham);
        TheBlock nextBlock(const stepData& data, rmMatrixXd& psiGround);
                                                     // performs each DMRG step
        obsMatrixX_t obsChangeBasis(const obsMatrixX_t& mat) const;
                       // changes basis during calculation of observables stage
        FinalSuperblock createHSuperFinal(const stepData& data,
                                          const rmMatrixXd& psiGround,
                                          int skips) const;
    
    private:
        Eigen::MatrixXd hS;                                // block Hamiltonian
        std::vector<Eigen::MatrixXd> rhoBasisH2;
                                     // density-matrix-basis coupling operators
        
        Eigen::MatrixXd changeBasis(const Eigen::MatrixXd& mat) const;
                   // represents operators in the basis of the new system block
};

#endif
