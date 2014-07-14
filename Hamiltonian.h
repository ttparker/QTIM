#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"
#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)

typedef std::vector<MatrixD_t, Eigen::aligned_allocator<MatrixD_t>> vecMatD_t;

class Hamiltonian
{
    public:
        int lSys;                                      // current system length
        
        Hamiltonian();
        void setParams(const std::vector<double>& couplingConstants, int lSys);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        std::vector<double> couplingConstants;
        vecMatD_t siteBasisH2;                 // site-basis coupling operators
        MatrixD_t westSideH1,
              // single-site Hamiltonian for initially left-hand half of system
                  eastSideH1;
             // single-site Hamiltonian for initially right-hand half of system
        
        MatrixX_t blockSiteJoin(const std::vector<MatrixX_t>& rhoBasisH2) const,
                                           // appends free site to system block
                  siteSiteJoin(int m, int compm) const;
                                           // joins the two free sites together
    
    friend class TheBlock;
};

#endif
