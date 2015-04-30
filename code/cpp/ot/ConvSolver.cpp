#include "ConvSolver.h"
#include "Timer.h"
#include <cassert>

namespace ot
{

double ConvSolver::computeKL(const VectorXd& w0, const VectorXd& w1) const
{
    VectorXd l0 = w0.array() * w0.array().log();
    VectorXd l1 = w1.array() * w1.array().log();

    VectorXd Kw0 = applyKernel(w0);
    VectorXd Kw1 = applyKernel(w1);

    double sum = 0.0;
    sum += dotProduct(l0, Kw1);
    sum += dotProduct(l1, Kw0);
    return sum;
}

double ConvSolver::computeMarginalEntropy(const VectorXd& q, double beta) const
{
    // - \int_M q(x)^beta*log(q(x)^beta) dx
    VectorXd l = beta * q.array().log();
    VectorXd e = l.array().exp();
    return - dotProduct(e, l);
}

double ConvSolver::computeDistance(VectorXd p0, VectorXd p1, int verbose) const
{
    VectorXd w0, w1;
    double val = computeDistance(p0, p1, w0, w1, verbose);
    return val;
}

double ConvSolver::computeDistance(VectorXd  p0, VectorXd  p1,
                                   VectorXd& w0, VectorXd& w1,
                                   int verbose) const
{
    // sanity check
    assert(mArea.size() > 0);

    assert(p0.size() == p1.size());
    assert(p0.size() == mArea.size());

    assert(std::abs((p0.array() * mArea.array()).sum() - 1.0) < 1e-12); 
    assert(std::abs((p1.array() * mArea.array()).sum() - 1.0) < 1e-12); 

    std::vector<double> timer;
    if (verbose > 0) Timer::start(timer, COLOR_BLUE, "distance");

    // init
    unsigned n = mArea.size();

    // smooth a bit to make less noisy
    p0.array() += mOpt.epsilon;
    p1.array() += mOpt.epsilon;

    w0 = VectorXd::Constant(n, 1.);
    w1 = VectorXd::Constant(n, 1.);
    
    VectorXd tmp;
    double dist = 0.0;
    unsigned iter = 0;
    for ( ; iter < mOpt.maxIters; ++iter)
    {
        // update w0
        tmp = applyKernel(w1);
        clamp(tmp, mOpt.epsilon);
        w0 = p1.array() / tmp.array();

        // update w1
        tmp = applyKernel(w0);
        clamp(tmp, mOpt.epsilon);
        w1 = p0.array() / tmp.array();

        // update dist
        double prev = dist;
        dist = computeKL(w0, w1);

        // stop criterion
        double res = std::abs(prev - dist);
        if (verbose > 1) std::cout << "Iter" << iter << ": " << res << std::endl;
        if (iter > 0 and res < mOpt.tolerance) break;
    }
    dist = std::sqrt(mOpt.gamma * std::max(dist, 0.0));

    if (verbose > 0) Timer::stop(timer, COLOR_BLUE);
    if (verbose > 1) 
        std::cout 
        << "dist = " << dist << " ; "
        << "converged in " 
        << iter << " / " << mOpt.maxIters 
        << " iterations" << std::endl;

    return dist;
}

void ConvSolver::computeBarycenter(std::vector<VectorXd> p, VectorXd alpha, VectorXd& q, bool useSharpening, int verbose) const
{
    // sanity check
    assert(p.size() > 0);
    assert((int)p.size() == alpha.size());

    assert(mArea.size() > 0);
    for (unsigned i = 0; i < p.size(); ++i) 
    {
        assert(p[i].size() == mArea.size());
        assert(std::abs((p[i].array() * mArea.array()).sum() - 1.0) < 1e-12); 
    }
    
    std::vector<double> timer;
    if (verbose > 0) Timer::start(timer, COLOR_BLUE, "barycenter");

    // init
    unsigned k = p.size();
    unsigned n = mArea.size();

    alpha /= alpha.sum();
    for (int i = 0; i < int(k); ++i) p[i].array() += mOpt.epsilon;

    q = VectorXd::Constant(n, 1.);
    std::vector<VectorXd> v(k, VectorXd::Constant(n, 1.));
    std::vector<VectorXd> w(k, VectorXd::Constant(n, 1.));
    
    unsigned iter = 0;
    for ( ; iter < mOpt.maxIters; ++iter)
    {
        // update w
        #pragma omp parallel for
        for (int i = 0; i < int(k); ++i)
        {
            VectorXd tmp = applyKernel(v[i]);
            clamp(tmp, mOpt.epsilon);
            w[i] = p[i].array() / tmp.array();
        }

        // compute auxiliar d
        std::vector<VectorXd> d(k);
        #pragma omp parallel for
        for (int i = 0; i < int(k); ++i)
        {
            d[i] = v[i].array() * applyKernel(w[i]).array();
            clamp(d[i], mOpt.epsilon);
        }

        // update barycenter q
        VectorXd prev = q;
        q = VectorXd::Zero(n);
        for (int i = 0; i < int(k); ++i) 
        {
            q = q.array() + alpha[i] * d[i].array().log();
        }
        q = q.array().exp();

        // Optional
        if (useSharpening and iter > 0) sharpen(q);

        // update v
        #pragma omp parallel for
        for (int i = 0; i < int(k); ++i) 
        {
            v[i] = v[i].array() * (q.array() / d[i].array());
        }

        // stop criterion
        VectorXd diff = prev - q;
        double res = dotProduct(diff, diff);
        if (verbose > 0) std::cout << green << "Iter " << white << iter << ": " << res << std::endl;
        if (iter > 1 and res < mOpt.tolerance) break;
    }

    if (verbose > 0) Timer::stop(timer, COLOR_BLUE);
    if (verbose > 1) 
        std::cout 
        << "Converged in " 
        << iter << " / " << mOpt.maxIters 
        << " iterations" << std::endl;
}

void ConvSolver::sharpen(VectorXd& q) const
{
    const double tol = mOpt.tolerance;

    // q has low entropy
    double val = computeMarginalEntropy(q, 1.0) - mOpt.upperEntropy;
    if (std::abs(val) < tol) return;

    // binary search for 'attempts' times
    double alpha0 = 0.8;
    double alpha1 = 1.2;
    unsigned attempt = 100;
    double beta = binarySearch(q, alpha0, alpha1, attempt, tol);

    if (attempt == 0) return; // binary search failed
    q = q.array().pow(beta);  // finally, sharpening
}

double ConvSolver::binarySearch(const VectorXd& q, double x0, double x1, unsigned& attempt, double tol) const
{
    if (attempt == 0) return 1.0; // abort
    attempt--;

    double f0 = computeMarginalEntropy(q, x0) - mOpt.upperEntropy;
    if (std::abs(f0) < tol) return x0; // x0 is root
    
    double f1 = computeMarginalEntropy(q, x1) - mOpt.upperEntropy;
    if (std::abs(f1) < tol) return x1; // x1 is root

    if (std::abs(f1-f0) < tol) { attempt = 0; return 1.0; } // abort

    double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    return binarySearch(q, x1, x2, attempt, tol);
}

} // namespace ot
