#ifndef CONV_OT_OPTIONS_H
#define CONV_OT_OPTIONS_H

namespace ot
{

class Options
{
public:
    double   upperEntropy; // entropy upper bound
    double   tolerance;    // convergence tolerance
    double   diffIters;    // number of diffusion iterations
    unsigned maxIters;     // number of projection iterations
    double   epsilon;      // numerical zero
    double   gamma;        // sharpening parameter

    Options()
    : upperEntropy(1.0)
    , tolerance(1e-6)
    , diffIters(10)
    , maxIters(1000)
    , epsilon(1e-20)
    , gamma(0.0)
    { }

    Options(const Options& rhs) 
    { 
        copy(rhs); 
    }

    Options& operator = (const Options& rhs) 
    { 
        copy(rhs); return *this; 
    }

    void copy(const Options& rhs)
    {
        upperEntropy = rhs.upperEntropy;
        tolerance = rhs.tolerance;
        diffIters = rhs.diffIters;
        maxIters = rhs.maxIters;
        epsilon = rhs.epsilon;
        gamma = rhs.gamma;
    }
};

} // namespace ot

#endif // CONV_OT_OPTIONS_H
