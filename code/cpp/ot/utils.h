#ifndef OT_UTILS_H
#define OT_UTILS_H

#include "types.h"
#include "LinearSolver.h"

namespace ot
{

inline
void clamp(VectorXd& q, double zero)
{
    for (int i = 0; i < int(q.size()); ++i) 
    {
        if (q[i] > zero or q[i] < -zero) continue;
        if (q[i] >= 0.0) q[i] =  zero;
        else             q[i] = -zero;
    }
}

inline
double dot_product(const VectorXd& area, const VectorXd& x, const VectorXd& y)
{
    VectorXd tmp = area.array() * y.array();
    return x.dot(tmp);
}

inline
VectorXd apply_kernel(const LinearSolver& lsolver, 
                      const VectorXd& area, 
                      const VectorXd& rhs,
                      unsigned iters = 1)
{
    VectorXd x = rhs;
    for (unsigned i = 0; i < iters; ++i)
    {
        VectorXd b = area.array() * x.array();
        lsolver.solvePosDef(b, x);
    }
    return x;
}

} // namespace ot

#endif // OT_UTILS_H
