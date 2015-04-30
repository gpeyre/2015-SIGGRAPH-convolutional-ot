#include "LinearSolver.h"

bool 
LinearSolver::factorizePosDef(const SparseMatrix& A)
{
    mPosDefSolver.compute(A);
    return (mPosDefSolver.info() == Eigen::Success);
}

bool 
LinearSolver::solvePosDef(const VectorXd& B, VectorXd& X) const
{
    X.resize(B.size());
    X = mPosDefSolver.solve(B);
    return (mPosDefSolver.info() == Eigen::Success);
}
