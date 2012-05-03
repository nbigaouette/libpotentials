#ifndef INC_VECTORS_hpp
#define INC_VECTORS_hpp

#include <cmath>

#ifdef YDEBUG
#include <stdexcept>
#include <limits>  // http://www.cplusplus.com/reference/std/limits/numeric_limits/
#endif // #ifdef YDEBUG

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
Double Get_Distance_Squared(const Double r1[3], const Double r2[3])
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
Double Get_Distance(const Double r1[3], const Double r2[3])
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
        const Double r1[3], const Double r2[3],
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

// **************************************************************
template <class Double>
void Rodrigues_Rotation(const Double v[3], const Double k[3], const Double theta,
                        Double vrot[3])
/**
 * Rotate in 3D vector "v" around axis "k" of angle "theta".
 * See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
 * @param   v       Vector to rotate
 * @param   k       Unit vector around which the rotation in to be taken
 * @param   theta   Angle of rotation
 * @param   vrot    Rotated vector to return
 */
{
    const Double sinTheta = std::sin(theta);
    const Double cosTheta = std::cos(theta);
    const Double k_dot_v  = Vector_Dot_Product(k, v);
    Double k_cross_v[3];
    Vector_Cross_Product(k_cross_v, k, v);
    const Double one_minus_cosTheta = Double(1.0) - cosTheta;

#ifdef YDEBUG
    // Vector to rotate must not be 0
    const Double v_norm = Vector_Length<Double>(v);
    if (v_norm < 100.0*std::numeric_limits<Double>::min())
    {
        throw std::runtime_error("Vector to rotate cannot be 0!");
    }

    // Rotation axis k MUST be unitary!
    const Double k_norm = Vector_Length<Double>(k);
    if (not (1.0-1.0e-7 < k_norm && k_norm < 1.0+1.0e-7))
    {
        throw std::runtime_error("Rotation axis MUST be unitary!");
    }
#endif // #ifdef YDEBUG

    for (int d = 0 ; d < 3 ; d++)
    {
        vrot[d] = (v[d] * cosTheta) + (k_cross_v[d] * sinTheta) + (k[d] * k_dot_v * one_minus_cosTheta);
    }
}

#endif // INC_VECTORS_hpp

// ********** End of file ***************************************
