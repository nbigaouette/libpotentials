#ifndef INC_VECTORS_hpp
#define INC_VECTORS_hpp

#include <cmath>

// **************************************************************
//                      Vector operations
// **************************************************************
template <class Double>
void Vector_Cross_Product(Double a[3], const Double b[3], const Double c[3])
{
    a[0] = b[1]*c[2] - b[2]*c[1];
    a[1] = b[2]*c[0] - b[0]*c[2];
    a[2] = b[0]*c[1] - b[1]*c[0];
}

template <class Double>
Double Vector_Dot_Product(const Double a[3], const Double b[3])
{
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

template <class Double>
Double Vector_Length_Squared(const Double a[3])
{
    return Vector_Dot_Product(a, a);
}

template <class Double>
Double Vector_Length(const Double a[3])
{
    return std::sqrt(Vector_Length_Squared(a));
}

template <class Double>
void Vector_Times_Scalar(Double a[3], const Double b[3], const Double c)
{
    a[0] = b[0] * c;
    a[1] = b[1] * c;
    a[2] = b[2] * c;
}

template <class Double>
void Vector_Add(Double a[3], const Double b[3], const Double c[3])
{
    a[0] = b[0] + c[0];
    a[1] = b[1] + c[1];
    a[2] = b[2] + c[2];
}

template <class Double>
void Vector_Substract(Double a[3], const Double b[3], const Double c[3])
{
    a[0] = b[0] - c[0];
    a[1] = b[1] - c[1];
    a[2] = b[2] - c[2];
}

template <class Double>
Double Vector_Substract_and_Length_Squared(Double a[3], const Double b[3])
{
    return (
              (a[0] - b[0])*(a[0] - b[0])
            + (a[1] - b[1])*(a[1] - b[1])
            + (a[2] - b[2])*(a[2] - b[2])
        );
}

// **************************************************************
template <class Double>
void Get_r21(const Double r1[3], const Double r2[3], Double r21[3])
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1   Input:  Position r1 [any units]
 * @param  r2   Input:  Position r2 [any units]
 * @param  r21  Output: Vector from to 2 to 1 [same units]
 */
{
    for (int d = 0 ; d < 3 ; d++) r21[d] = r1[d] - r2[d];
}

// **************************************************************
template <class Double>
Double Get_Distance_Squared(Double r1[3], Double r2[3])
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1    Position r1 [any units]
 * @param  r2    Position r2 [any units]
 * @return r212  Length of vector from to 2 to 1 squared [same units^2]
 */
{
    Double r212 = 0.0;
    Double r21[3];
    Get_r21(r1, r2, r21);
    for (int d = 0 ; d < 3 ; d++) r212 += r21[d]*r21[d];
    return r212;
}

// **************************************************************
template <class Double>
Double Get_Distance(Double r1[3], Double r2[3])
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1    Position r1 [any units]
 * @param  r2    Position r2 [any units]
 * @return r212  Length of vector from to 2 to 1 [same units]
 */
{
    return std::sqrt(Get_Distance_Squared(r1, r2));
}

// **************************************************************
template <class Double>
void set_vector_between_particles(
        Double r1[3], Double r2[3],
        Double r21[3], Double &r212, Double &r, Double &one_over_r)
/**
 * r21 is the vector groing from position r2 to position r1
 * r21 = r1 - r2
 * r2 + r21 = r1
 *
 *     2
 *    / \
 *   /   \ dr
 *  /    _\/
 * /------>1
 *
 * @param  r1   Input:  Position 1 [any units]
 * @param  r2   Input:  Position 2 [any units]
 * @param  r21  Output: Vector from to 2 to 1 [same units]
 * @param  r212 Output: Length squared of vector r12 [same units^2]
 * @param  r    Output: Length of vector r12 [same units]
 * @param  one_over_r Output: Inverse of the length of vector r12 [same units^-1]
 */
{
    Get_r21(r1, r2, r21);
    r212       = Get_Distance_Squared(r1, r2);
    r          = std::sqrt(r212);
    one_over_r = Double(1.0) / r;
}

#endif // INC_VECTORS_hpp

// ********** End of file ***************************************
