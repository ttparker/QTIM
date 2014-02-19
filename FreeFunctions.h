#ifndef FREEFUNCTIONS_H
#define FREEFUNCTIONS_H

#include "EffectiveHamiltonian.h"

void oneSiteExpValues(const MatrixDd& oneSiteOp, int rangeOfObservables,
					  int currentLSys, EffectiveHamiltonian& hSuperFinal,
					  std::vector<TheBlock>& blocks, std::ofstream& fileout);
void twoSiteExpValues(const MatrixDd& firstTwoSiteOp,
					  const MatrixDd& secondTwoSiteOp, int rangeOfObservables,
					  int currentLSys, EffectiveHamiltonian& hSuperFinal,
					  std::vector<TheBlock>& blocks, std::ofstream& fileout);
void modifyHamParams(int trial = 0);

#endif
