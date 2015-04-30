#ifndef OT_CONV_SOLVER_H
#define OT_CONV_SOLVER_H

#include "utils.h"
#include "Options.h"

namespace ot
{

class ConvSolver
{
public:
    ConvSolver(const Options& opt, const VectorXd& area, const LinearSolver& solver)
    : mOpt(opt)
    , mArea(area)
    , mSolver(solver)
    { 
    }

    double computeKL(const VectorXd& w0, const VectorXd& w1) const;
    double computeMarginalEntropy(const VectorXd& q, double beta = 1.0) const;

    double computeDistance(VectorXd p0, VectorXd p1, int verbose = 0) const;
    double computeDistance(VectorXd p0, VectorXd p1, VectorXd& w0, VectorXd& w1, int verbose = 0) const;

    void computeBarycenter(std::vector<VectorXd> p, VectorXd alpha, VectorXd& q, bool useSharpening, int verbose = 0) const;

    void   sharpen(VectorXd& q) const;
    double binarySearch(const VectorXd& q, double x0, double x1, unsigned& attempt, double tol) const;

    // utils wrappers //

    inline double dotProduct(const VectorXd& x, const VectorXd& y) const
    { return dot_product(mArea, x, y); }

    inline VectorXd applyKernel(const VectorXd& rhs) const 
    { return apply_kernel(mSolver, mArea, rhs, mOpt.diffIters); }

    inline       Options& options()       { return mOpt; }
    inline const Options& options() const { return mOpt; }

private:
    // members
    Options mOpt;
    const VectorXd& mArea;
    const LinearSolver& mSolver;
};

} // namespace ot

#endif // OT_CONV_SOLVER_H
