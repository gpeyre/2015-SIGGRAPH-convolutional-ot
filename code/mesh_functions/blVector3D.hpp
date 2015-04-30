#ifndef BL_VECTOR3D_HPP
#define BL_VECTOR3D_HPP

#include <cmath>
#include <iostream>
#include <cstdlib>

//-------------------------------------------------------------------
// CLASS:           blVector3d
// BASE CLASS:      None
// PURPOSE:         A simple 3d vector class
// AUTHOR:          Vincenzo Barbato
//                  http://www.barbatolabs.com
//                  navyenzo@gmail.com
// LISENSE:         MIT-LICENCE
//                  http://www.opensource.org/licenses/mit-license.php
// DEPENDENCIES:    std::abs -- To calculate norms
//                  std::sqrt -- To calculate magnitude
//                  iostream -- To output to the console
// NOTES:           - This 3d vector class has an additional operator
//                    defined in this file that allows multiplying a
//                    scalar by a 3d vector
//                  - CrossProduct -- Defined in this file to take the
//                                    cross product of two vectors
// DATE CREATED:    Mar/30/2010
// DATE UPDATED:
//-------------------------------------------------------------------


//-------------------------------------------------------------------
// Includes and libs needed for this file
//-------------------------------------------------------------------
//-------------------------------------------------------------------


//-------------------------------------------------------------------
// Enums used for this file and sub-files
//-------------------------------------------------------------------
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
class blVector3d
{
public: // Constructors and destructors

    // Default constructor
    blVector3d(const blType& X = 0,const blType& Y = 0,const blType& Z = 0);

    // Constructor from a Vector of different type
    template<typename blType2>
    blVector3d(const blVector3d<blType2> &Vector);

    // Destructor
    ~blVector3d();

public: // Overloaded operators

    const bool                  operator==(const blVector3d<blType>& Vector)const;
    const bool                  operator!=(const blVector3d<blType>& Vector)const;

    const blVector3d<blType>    operator-()const;
    blVector3d<blType>&         operator+=(const blVector3d<blType>& Vector);
    blVector3d<blType>&         operator-=(const blVector3d<blType>& Vector);
    blVector3d<blType>&         operator*=(const blType& Scalar);
    blVector3d<blType>&         operator/=(const blType& Scalar);

    const blVector3d<blType>    operator+(const blVector3d<blType>& Vector)const;
    const blVector3d<blType>    operator-(const blVector3d<blType>& Vector)const;
    const blVector3d<blType>    operator*(const blType& Scalar)const;
    const blVector3d<blType>    operator/(const blType& Scalar)const;
    
    blType &operator[](int i) {
        if (i == 0) return x;
        else if (i == 1) return y;
        else if (i == 2) return z;
        else {
            std::cout << "Bad vector component.\n";
            exit(0);
        }
    }

    // Dot product
    const blType                operator*(const blVector3d<blType>& Vector)const;

public: // Public functions

    // Vector magnitude
    const blType            GetMagnitude()const;

    // Functions to get a normilized
    // Vector and to get a Vector
    // perpendicular to this one
    const blVector3d<blType> GetPerpendicularUnitVector()const;
    const blVector3d<blType> GetNormalizedVector()const;

    // Function used to get a vector
    // with inverted elements
    const blVector3d<blType> GetVectorWithInvertedElements()const;

    // Function used to normalize this Vector
    void                    NormalizeVector();

    // Functions to calculate norms
    const blType            Norm1()const;
    const blType            Norm2()const;
    const blType            NormInf()const;

public: // Public variables

    // The Vector's x,y and z coordinates
    blType                  x;
    blType                  y;
    blType                  z;
};
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType>::blVector3d(const blType& X,const blType& Y,const blType& Z)
{
    x = X;
    y = Y;
    z = Z;
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
template<typename blType2>
inline blVector3d<blType>::blVector3d(const blVector3d<blType2>& Vector)
{
    x = (blType)Vector.x;
    y = (blType)Vector.y;
    z = (blType)Vector.z;
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType>::~blVector3d(void)
{
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::GetVectorWithInvertedElements()const
{
    return blVector3d<blType>(blType(1)/x,blType(1)/y,blType(1)/z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const bool blVector3d<blType>::operator==(const blVector3d<blType>& Vector)const
{
    // Check the components
    if(x != Vector.x)
        return false;
    if(y != Vector.y)
        return false;
    if(z != Vector.z)
        return false;

    // By now we know the components
    // match so the vectors are equal
    return true;
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const bool blVector3d<blType>::operator!=(const blVector3d<blType>& Vector)const
{
    // Check the components
    if(x != Vector.x)
        return true;
    if(y != Vector.y)
        return true;
    if(z != Vector.z)
        return true;

    // By now we know the components
    // match so the vectors are equal
    return false;
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::operator-()const
{
    return blVector3d<blType>(-x,-y,-z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType>& blVector3d<blType>::operator+=(const blVector3d<blType> &Vector)
{
    x += Vector.x;
    y += Vector.y;
    z += Vector.z;

    return (*this);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType>& blVector3d<blType>::operator-=(const blVector3d<blType> &Vector)
{
    x -= Vector.x;
    y -= Vector.y;
    z -= Vector.z;

    return (*this);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType>& blVector3d<blType>::operator*=(const blType& Scalar)
{
    x *= Scalar;
    y *= Scalar;
    z *= Scalar;

    return (*this);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType>& blVector3d<blType>::operator/=(const blType& Scalar)
{
    x /= Scalar;
    y /= Scalar;
    z /= Scalar;

    return (*this);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType> operator*(const blType& Scalar,
                                    const blVector3d<blType>& Vector)
{
    return blVector3d<blType>(Scalar*Vector.x,Scalar*Vector.y,Scalar*Vector.z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::operator*(const blType& Scalar)const
{
    return blVector3d<blType>(x * Scalar,y * Scalar,z * Scalar);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::operator/(const blType& Scalar)const
{
    return blVector3d<blType>(x / Scalar,y / Scalar, z / Scalar);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::operator+(const blVector3d<blType>& Vector)const
{
    return blVector3d<blType>(x + Vector.x,y + Vector.y, z + Vector.z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::operator-(const blVector3d<blType>& Vector)const
{
    return blVector3d<blType>(x - Vector.x,y - Vector.y,z - Vector.z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blType blVector3d<blType>::operator*(const blVector3d<blType>& Vector)const
{
    return (x * Vector.x + y * Vector.y + z * Vector.z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blType blVector3d<blType>::GetMagnitude()const
{
    return std::sqrt(blType(x*x + y*y + z*z));
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::GetPerpendicularUnitVector()const
{
    // Find a component of this Vector that is non zero
    // and if unsuccessfull, just return the zero Vector
    //
    // P1*P2 = x1*x2 + y1*y2 + z1*z2 = 0

    if(x != 0)
    {
        if(y == 0)
            return blVector3d<blType>(0,1,0);
        else
        {
            if(z == 0)
            {
                return blVector3d<blType>(0,0,1);
            }
            else
            {
                blVector3d<blType> Vector(1,1,(x + y)/z);
                Vector.NormalizeVector();

                return Vector;
            }
        }
    }
    else if(y != 0 || z != 0)
    {
        return blVector3d<blType>(1,0,0);
    }

    // The Vector was a zero Vector, therefore any Vector would be perpendicular,
    // thus we are just going to return the zero Vector
    // Error
    return blVector3d<blType>(0,0,0);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blVector3d<blType> blVector3d<blType>::GetNormalizedVector()const
{
    blType Magnitude = this->GetMagnitude();
    return blVector3d<blType>(x/Magnitude,y/Magnitude,z/Magnitude);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline void blVector3d<blType>::NormalizeVector()
{
    blType Magnitude = GetMagnitude();

    x /= Magnitude;
    y /= Magnitude;
    z /= Magnitude;
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blType blVector3d<blType>::Norm1()const
{
    return (std::abs(x) + std::abs(y) + std::abs(z));
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blType blVector3d<blType>::Norm2()const
{
    return std::sqrt(x*x + y*y + z*z);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline const blType blVector3d<blType>::NormInf()const
{
    if(x > y)
    {
        if(x > z)
            return x;
        else
            return z;
    }
    else
    {
        if(y > z)
            return y;
        else
            return z;
    }
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline blVector3d<blType> CrossProduct(const blVector3d<blType> V1,
                                       const blVector3d<blType> V2)
{
    return blVector3d<blType>(V1.y*V2.z - V1.z*V2.y,
                              -V1.x*V2.z + V1.z*V2.x,
                              V1.x*V2.y - V1.y*V2.x);
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template<typename blType>
inline std::ostream& operator<<(std::ostream& os,const blVector3d<blType>& Vector)
{
    os << "(" << Vector.x << "," << Vector.y << "," << Vector.z << ")";
    return os;
}
//-------------------------------------------------------------------


#endif // BL_VECTOR3D_HPP
