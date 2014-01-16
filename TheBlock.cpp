#include "d.h"
#include "main.h"
#include "Hamiltonian.h"
#include "TheBlock.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const MatrixXd& hS,
				   const std::vector<MatrixXd>& rhoBasisH2)
				   : hS(hS), rhoBasisH2(rhoBasisH2), m(m) {};

TheBlock::TheBlock(const Hamiltonian& ham, int mMaxIn) : hS(ham.h1), m(d)
{
	mMax = mMaxIn;
	rhoBasisH2.assign(ham.h2.begin(),
					  ham.h2.begin() + ham.couplingConstants.size());
};

TheBlock TheBlock::nextBlock(const Hamiltonian& ham, bool infiniteStage,
							 TheBlock& compBlock)
{
    MatrixXd hSprime = kp(hS, Id_d)
                       + ham.blockSiteJoin(rhoBasisH2)
                       + kp(Id(m), ham.h1);		// expanded system block
	std::vector<MatrixXd> tempRhoBasisH2;
	int indepCouplingOperators = ham.couplingConstants.size();
	tempRhoBasisH2.reserve(indepCouplingOperators);
	int md = m * d;
	if(md <= mMax)
	{ // if near edge of system, no truncation necessary so skip DMRG algorithm
		for(auto op = ham.h2.begin(), end = ham.h2.begin() + indepCouplingOperators;
			op != end; op++)
			tempRhoBasisH2.push_back(kp(Id(m), *op));
		return TheBlock(md, hSprime, tempRhoBasisH2);
	};
	SelfAdjointEigenSolver<MatrixXd> hSuperSolver(infiniteStage ?
												// find superblock eigenstates
												  MatrixXd(kp(hSprime, Id(md))
												  + ham.siteSiteJoin(m, m)
												  + kp(Id(md), hSprime)) :
												  MatrixXd(kp(hSprime, Id(compBlock.m * d))
								  			      + ham.siteSiteJoin(m, compBlock.m)
												  + kp(Id(md * compBlock.m), ham.h1)
											      + kp(Id(md), ham.blockSiteJoin(compBlock.rhoBasisH2))
												  + kp(kp(Id(md), compBlock.hS), Id_d)));
    rmMatrixXd psiGround = hSuperSolver.eigenvectors().col(0);	// ground state
    psiGround.resize(md, (infiniteStage ? m : compBlock.m) * d);
    SelfAdjointEigenSolver<MatrixXd> rhoSolver(psiGround * psiGround.adjoint());
											// find density matrix eigenstates
	primeToRhoBasis = rhoSolver.eigenvectors().rightCols(mMax);
											// construct change-of-basis matrix
	for(auto op = ham.h2.begin(),
			 end = ham.h2.begin() + indepCouplingOperators; op != end; op++)
		tempRhoBasisH2.push_back(changeBasis(kp(Id(m), *op)));
	return TheBlock(mMax, changeBasis(hSprime), tempRhoBasisH2);
								// save expanded-block operators in new basis
};

std::pair<Eigen::MatrixXd, int> TheBlock::createHSuperFinal(const Hamiltonian& ham)
{
    return std::make_pair(MatrixXd(kp(hS, Id(d * m * d))
								   + kp(ham.blockSiteJoin(rhoBasisH2), Id(m * d))
								   + kp(kp(Id(m), ham.h1), Id(m * d))
								   + ham.siteSiteJoin(m, m)
								   + kp(Id(m * d * m), ham.h1)
								   + kp(Id(m * d), ham.blockSiteJoin(rhoBasisH2))
								   + kp(kp(Id(m * d), hS), Id_d)),
								   m);
};

MatrixXd TheBlock::changeBasis(const MatrixXd& mat)
{
	return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
