#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "types.h"

#include <Eigen/SparseCholesky>

class LinearSolver
{
public:
    bool factorizePosDef(const SparseMatrix& A);

    bool solvePosDef(const VectorXd& B, VectorXd& X) const;

private:
    Eigen::SimplicialLDLT<SparseMatrix> mPosDefSolver;
};


#endif // LINEAR_SOLVER_H
