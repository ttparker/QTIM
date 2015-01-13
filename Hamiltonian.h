#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"

struct effectiveHams       // each site's effective block and bond Hamiltonians
{
    int m;                                  // number of states stored in block
    rmMatrixX_t hS;                                        // block Hamiltonian
    std::vector<MatrixX_t> rhoBasisH2; // density-matrix-basis coupling operators
};

class TheBlock;

class Hamiltonian
{
    public:
        Hamiltonian();
        void setParams(const std::vector<double>& couplingConstants);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        std::vector<double> couplingConstants;
        std::vector<MatrixD_t, Eigen::aligned_allocator<MatrixD_t>> siteBasisH2;
                                               // site-basis coupling operators
        MatrixD_t h1;                                // single-site Hamiltonian
        
        rmMatrixX_t EDNNCoupling(const std::vector<MatrixX_t>& rhoBasisH2) const;
        Eigen::Matrix<hamScalar, Eigen::Dynamic, 1>
            act(const effectiveHams& blockParts,
                const effectiveHams& compBlockParts, rmMatrixX_t x) const;
        rmMatrixX_t NNCoupling(TheBlock * block) const;
                                    // specify the block-adjacent-free-site
                                    // couplings to be projected into new basis
    
    friend class TheBlock;
};

#endif
