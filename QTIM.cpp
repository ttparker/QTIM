#include "Hamiltonian.h"

#define j couplingConstants[0]
#define h oneSiteConstants[0]
#define sigmax h2[0]
#define rhoBasisSigmax rhoBasisH2[0]

using namespace Eigen;

Hamiltonian::Hamiltonian(int lSys,
						 const std::vector<double>& couplingConstants,
						 const std::vector<double>& oneSiteConstants)
						 : lSys(lSys), couplingConstants(couplingConstants)
{
	h2.resize(couplingConstants.size());
	sigmax << 0., 1.,
			  1., 0.;								 // define Pauli matrices
	h1 << -h, 0.,
		  0., h;
};

MatrixXd Hamiltonian::blockSiteJoin(const std::vector<MatrixXd>& rhoBasisH2) const
{
	return j * kp(rhoBasisSigmax, sigmax);
};

MatrixXd Hamiltonian::siteSiteJoin(int m1, int m2) const
{
	return j * kp(kp(Id(m1), sigmax), kp(Id(m2), sigmax));
};
