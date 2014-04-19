#ifndef FREEFUNCTIONS_H
#define FREEFUNCTIONS_H

#include "EffectiveHamiltonian.h"

rmMatrixXd randomSeed(const TheBlock& leftBlock, const TheBlock& rightBlock);
                                            // outputs random normalized vector
void reflectPredictedPsi(rmMatrixXd& psiGround, const TheBlock& bigBlock,
                         const TheBlock& littleBlock);
                                               // when you reach edge of system
Eigen::VectorXd oneSiteExpValues(const MatrixDd& oneSiteOp,
                                 int rangeOfObservables, int lSys,
                                 EffectiveHamiltonian& hSuperFinal,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks,
                                 std::ofstream& fileout);
Eigen::MatrixXd twoSiteExpValues(const MatrixDd& firstTwoSiteOp,
                                 const MatrixDd& secondTwoSiteOp,
                                 int rangeOfObservables, int lSys,
                                 EffectiveHamiltonian& hSuperFinal,
                                 std::vector<TheBlock>& leftBlocks,
                                 std::vector<TheBlock>& rightBlocks,
                                 std::ofstream& fileout);

#endif
