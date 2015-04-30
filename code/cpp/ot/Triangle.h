#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "types.h"
#include <cmath>

class Triangle
{
public:
    static double area(const Vec3& a, const Vec3& b, const Vec3& c);
    // area based on points

    static double angle(const Vec3& a, const Vec3& b, const Vec3& c);
    // angle at a based on points

    static double cotan(const Vec3& a, const Vec3& b, const Vec3& c);
    // cotan at a based on points
};

//
// inline definitions
//

inline double 
Triangle::area(const Vec3& a, const Vec3& b, const Vec3& c) 
{
    return 0.5 * (b-a).cross(c-b).norm();
}

inline double 
Triangle::angle(const Vec3& a, const Vec3& b, const Vec3& c)
{
    Vec3 p01 = b - a;
    Vec3 p02 = c - a;
    return std::atan2(p01.cross(p02).norm(), p01.dot(p02));
}

inline double 
Triangle::cotan(const Vec3& a, const Vec3& b, const Vec3& c)
{
    Vec3 p01 = b - a;
    Vec3 p02 = c - a;
    return p01.dot(p02) / p01.cross(p02).norm();    
}

#endif // TRIANGLE_H
