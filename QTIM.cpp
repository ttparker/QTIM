#include "Hamiltonian.h"

#define j couplingConstants[0]
#define h oneSiteConstants[0]
#define sigmax h2[0]
#define rhoBasisSigmax rhoBasisH2[0]

using namespace Eigen;

Hamiltonian::Hamiltonian()
{
	h2.resize(1);
	sigmax << 0., 1.,
			  1., 0.;								 // define Pauli matrices
};

void Hamiltonian::setParams(int lSysIn, const std::vector<double>& couplingConstantsIn, double hIn)
{
    lSys = lSysIn;
    couplingConstants = couplingConstantsIn;
    h1 << -hIn, 0.,
          0.,   hIn;
}

MatrixXd Hamiltonian::blockSiteJoin(const std::vector<MatrixXd>& rhoBasisH2) const
{
	return j * kp(rhoBasisSigmax, sigmax);
};

MatrixXd Hamiltonian::siteSiteJoin(int m1, int m2) const
{
	return j * kp(kp(Id(m1), sigmax), kp(Id(m2), sigmax));
};
