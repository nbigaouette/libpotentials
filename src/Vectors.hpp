#ifndef INC_VECTORS_hpp
#define INC_VECTORS_hpp

#include <cmath>

// **************************************************************
//                      Vector operations
// **************************************************************
static inline void Vector_Cross_Product(double a[3], const double b[3], const double c[3])
{
    a[0] = b[1]*c[2] - b[2]*c[1];
    a[1] = b[2]*c[0] - b[0]*c[2];
    a[2] = b[0]*c[1] - b[1]*c[0];
}

static inline double Vector_Dot_Product(const double a[3], const double b[3])
{
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

static inline double Vector_Length_Squared(const double a[3])
{
    return Vector_Dot_Product(a, a);
}

static inline double Vector_Length(const double a[3])
{
    return sqrt(Vector_Length_Squared(a));
}

static inline void Vector_Times_Scalar(double a[3], const double b[3], const double c)
{
    a[0] = b[0] * c;
    a[1] = b[1] * c;
    a[2] = b[2] * c;
}

static inline void Vector_Add(double a[3], const double b[3], const double c[3])
{
    a[0] = b[0] + c[0];
    a[1] = b[1] + c[1];
    a[2] = b[2] + c[2];
}

static inline void Vector_Substract(double a[3], const double b[3], const double c[3])
{
    a[0] = b[0] - c[0];
    a[1] = b[1] - c[1];
    a[2] = b[2] - c[2];
}

static inline double Vector_Substract_and_Length_Squared(double a[3], const double b[3])
{
    return (
              (a[0] - b[0])*(a[0] - b[0])
            + (a[1] - b[1])*(a[1] - b[1])
            + (a[2] - b[2])*(a[2] - b[2])
        );
}

#endif // INC_VECTORS_hpp

// ********** End of file ***************************************
