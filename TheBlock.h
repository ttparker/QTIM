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
						   TheBlock& compBlock);	// performs each DMRG step
			// - the third argument is the environment block in the fDMRG stage
		std::pair<Eigen::MatrixXd, int> createHSuperFinal(const Hamiltonian& ham);

	private:
		Eigen::MatrixXd hS;								// block Hamiltonian
		std::vector<Eigen::MatrixXd> rhoBasisH2;
									// density-matrix-basis coupling operators
		int m;								// number of states stored in block
		static int mMax;				// max size of effective Hamiltonian
		Eigen::MatrixXd primeToRhoBasis;			// change-of-basis matrix

		Eigen::MatrixXd changeBasis(const Eigen::MatrixXd& mat);
				// represents operators in the basis of the new system block

	friend class EffectiveHamiltonian;
};
