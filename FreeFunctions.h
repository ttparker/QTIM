void halfSweep(std::vector<TheBlock>& blocks, int start,
			   const Hamiltonian& ham, bool infiniteStage);
												// perform half of a DMRG sweep
void oneSiteExpValues(const MatrixDd& oneSiteOp, int rangeOfObservables,
					  int currentLSys, EffectiveHamiltonian& hSuperFinal,
					  std::vector<TheBlock>& blocks, std::ofstream& fileout);
void twoSiteExpValues(const MatrixDd& firstTwoSiteOp,
					  const MatrixDd& secondTwoSiteOp, int rangeOfObservables,
					  int currentLSys, EffectiveHamiltonian& hSuperFinal,
					  std::vector<TheBlock>& blocks, std::ofstream& fileout);
