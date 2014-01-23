typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    rmMatrixXd;

class TheBlock
{
	public:
		TheBlock(int m = 0,
				 const Eigen::MatrixXd& hS = Eigen::MatrixXd(),
				 const std::vector<Eigen::MatrixXd>& rhoBasisH2
					= std::vector<Eigen::MatrixXd>());
		TheBlock(const Hamiltonian& ham, int mMax);
		TheBlock nextBlock(const Hamiltonian& ham, bool infiniteStage,
						   const TheBlock& compBlock); // performs each DMRG step
		std::pair<Eigen::MatrixXd, int> createHSuperFinal(const Hamiltonian& ham)
            const;

	private:
		Eigen::MatrixXd hS;								// block Hamiltonian
		std::vector<Eigen::MatrixXd> rhoBasisH2;
									// density-matrix-basis coupling operators
		int m;								// number of states stored in block
		static int mMax;				// max size of effective Hamiltonian
        static double lancTolerance;
        Eigen::MatrixXd primeToRhoBasis;			// change-of-basis matrix

		Eigen::MatrixXd changeBasis(const Eigen::MatrixXd& mat) const;
				// represents operators in the basis of the new system block

    friend void halfSweep(std::vector<TheBlock>& blocks, int start,
                          const Hamiltonian& ham, bool infiniteStage,
                          double lancTolerance);
	friend class EffectiveHamiltonian;
};
