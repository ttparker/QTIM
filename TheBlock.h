#ifndef THEBLOCK_H
#define THEBLOCK_H

#include "Hamiltonian.h"
#define Id_d MatrixDd::Identity()   // one-site identity matrix

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    rmMatrixXd;

class EffectiveHamiltonian;

class TheBlock
{
	public:
        static double lancTolerance;
        Eigen::MatrixXd primeToRhoBasis;            // change-of-basis matrix
        
		TheBlock(int m = 0,
				 const Eigen::MatrixXd& hS = Eigen::MatrixXd(),
				 const std::vector<Eigen::MatrixXd>& rhoBasisH2
					= std::vector<Eigen::MatrixXd>());
		TheBlock(const Hamiltonian& ham, int mMax);
		TheBlock nextBlock(const Hamiltonian& ham, bool exactDiag = true,
                           bool infiniteStage = true,
                           const TheBlock& compBlock = TheBlock(),
                           const TheBlock& beforeCompBlock = TheBlock());
                                                     // performs each DMRG step
        void randomSeed();                           // for iDMRG case
        void reflectPredictedPsi();            // when you reach edge of system
		EffectiveHamiltonian createHSuperFinal(const Hamiltonian& ham,
                                               const TheBlock& compBlock,
                                               int skips) const;
    
	private:
		Eigen::MatrixXd hS;								// block Hamiltonian
		std::vector<Eigen::MatrixXd> rhoBasisH2;
									// density-matrix-basis coupling operators
		int m;								// number of states stored in block
		static rmMatrixXd psiGround;
        static bool firstfDMRGStep;     // slight abuse of nomeclature
                                        //    - true during iDMRG as well
		static int mMax;				// max size of effective Hamiltonian
        
		Eigen::MatrixXd changeBasis(const Eigen::MatrixXd& mat) const;
				// represents operators in the basis of the new system block
    
	friend class EffectiveHamiltonian;
};

#endif
