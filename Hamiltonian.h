class Hamiltonian
{
	public:
		int lSys;									// current system length

		Hamiltonian(int lSys, const std::vector<double>& couplingConstants,
					const std::vector<double>& oneSiteConstants);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
	private:
		std::vector<double> couplingConstants;
		std::vector<MatrixDd, Eigen::aligned_allocator<MatrixDd>> h2;
                                            // site-basis coupling operators
		MatrixDd h1;								// single-site Hamiltonian

		Eigen::MatrixXd
			blockSiteJoin(const std::vector<Eigen::MatrixXd>& rhoBasisH2) const,
										// appends free site to system block
			siteSiteJoin(int m1, int m2) const;
										// joins the two free sites together

	friend class TheBlock;
    friend void modifyHamParams(int trial);
};
