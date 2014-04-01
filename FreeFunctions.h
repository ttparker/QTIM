#ifndef FREEFUNCTIONS_H
#define FREEFUNCTIONS_H

#include "EffectiveHamiltonian.h"

Eigen::VectorXd oneSiteExpValues(const MatrixDd& oneSiteOp,
                                 int rangeOfObservables, int currentLSys,
                                 EffectiveHamiltonian& hSuperFinal,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks,
                                 std::ofstream& fileout);
Eigen::MatrixXd twoSiteExpValues(const MatrixDd& firstTwoSiteOp,
                                 const MatrixDd& secondTwoSiteOp,
                                 int rangeOfObservables, int currentLSys,
                                 EffectiveHamiltonian& hSuperFinal,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks,
                                 std::ofstream& fileout);

#endif
