#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <limits>

#include "Constants.hpp"

#include "Vectors.hpp"

//const bool be_verbose = true;
const bool be_verbose = false;

using namespace libpotentials;


BOOST_AUTO_TEST_SUITE(vectors)

BOOST_AUTO_TEST_CASE(Test_Rodrigues_Rotation)
/**
 * Test function to rotate a vector around an axis using Rodrigues_Rotation()
 */
{
    {
        // Rotate (x,0,0) 90 degrees around z-axis == (0,x,0)
        const double to_rotate[3]       = {2.0, 0.0, 0.0};
        const double rotation_axis[3]   = {0.0, 0.0, 1.0};  // z-axis
        const double rotation_angle     = 90.0 * deg_to_rad;
        double result[3];
        const double answer[3]          = {0.0, 2.0, 0.0};
        Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result);

        //std::cout << "to_rotate = (" << to_rotate[0] << ", " << to_rotate[1] << ", " << to_rotate[2] << ")\n";
        //std::cout << "answer = (" << answer[0] << ", " << answer[1] << ", " << answer[2] << ")\n";
        //std::cout << "result = (" << result[0] << ", " << result[1] << ", " << result[2] << ")\n";

        BOOST_CHECK_SMALL(result[0] - answer[0], 1.0e-7);
        BOOST_CHECK_SMALL(result[1] - answer[1], 1.0e-7);
        BOOST_CHECK_SMALL(result[2] - answer[2], 1.0e-7);

        BOOST_CHECK_CLOSE(2.0,                              Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length(to_rotate),         Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length_Squared(to_rotate), Vector_Length_Squared(result),  0.0001);
    }

    {
        // Rotate (x,0,0) 45 degrees around z-axis == sqrt(2)/2 * (x,x,0)
        const double to_rotate[3]       = {2.0, 0.0, 0.0};
        const double rotation_axis[3]   = {0.0, 0.0, 1.0};  // z-axis
        const double rotation_angle     = 45.0 * deg_to_rad;
        double result[3];
        const double answer[3]          = {0.5*sqrt_2*to_rotate[0], 0.5*sqrt_2*to_rotate[0], 0.0};
        Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result);

        //std::cout << "to_rotate = (" << to_rotate[0] << ", " << to_rotate[1] << ", " << to_rotate[2] << ")\n";
        //std::cout << "answer = (" << answer[0] << ", " << answer[1] << ", " << answer[2] << ")\n";
        //std::cout << "result = (" << result[0] << ", " << result[1] << ", " << result[2] << ")\n";

        BOOST_CHECK_SMALL(result[0] - answer[0], 1.0e-7);
        BOOST_CHECK_SMALL(result[1] - answer[1], 1.0e-7);
        BOOST_CHECK_SMALL(result[2] - answer[2], 1.0e-7);

        BOOST_CHECK_CLOSE(2.0,                              Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length(to_rotate),         Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length_Squared(to_rotate), Vector_Length_Squared(result),  0.0001);
    }

    {
        // Rotate (x,0,0) -45 degrees around z-axis == sqrt(2)/2 * (x,-x,0)
        const double to_rotate[3]       = {2.0, 0.0, 0.0};
        const double rotation_axis[3]   = {0.0, 0.0, 1.0};  // z-axis
        const double rotation_angle     = -45.0 * deg_to_rad;
        double result[3];
        const double answer[3]          = {0.5*sqrt_2*to_rotate[0], -0.5*sqrt_2*to_rotate[0], 0.0};
        Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result);

        //std::cout << "to_rotate = (" << to_rotate[0] << ", " << to_rotate[1] << ", " << to_rotate[2] << ")\n";
        //std::cout << "answer = (" << answer[0] << ", " << answer[1] << ", " << answer[2] << ")\n";
        //std::cout << "result = (" << result[0] << ", " << result[1] << ", " << result[2] << ")\n";

        BOOST_CHECK_SMALL(result[0] - answer[0], 1.0e-7);
        BOOST_CHECK_SMALL(result[1] - answer[1], 1.0e-7);
        BOOST_CHECK_SMALL(result[2] - answer[2], 1.0e-7);

        BOOST_CHECK_CLOSE(2.0,                              Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length(to_rotate),         Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length_Squared(to_rotate), Vector_Length_Squared(result),  0.0001);
    }

    {
        // Rotate (0,y,0) 90 degrees around z-axis == (-y,0,0)
        const double to_rotate[3]       = {0.0, 2.0, 0.0};
        const double rotation_axis[3]   = {0.0, 0.0, 1.0};  // z-axis
        const double rotation_angle     = 90.0 * deg_to_rad;
        double result[3];
        const double answer[3]          = {-to_rotate[1], 0.0, 0.0};
        Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result);

        //std::cout << "to_rotate = (" << to_rotate[0] << ", " << to_rotate[1] << ", " << to_rotate[2] << ")\n";
        //std::cout << "answer = (" << answer[0] << ", " << answer[1] << ", " << answer[2] << ")\n";
        //std::cout << "result = (" << result[0] << ", " << result[1] << ", " << result[2] << ")\n";

        BOOST_CHECK_SMALL(result[0] - answer[0], 1.0e-7);
        BOOST_CHECK_SMALL(result[1] - answer[1], 1.0e-7);
        BOOST_CHECK_SMALL(result[2] - answer[2], 1.0e-7);

        BOOST_CHECK_CLOSE(2.0,                              Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length(to_rotate),         Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length_Squared(to_rotate), Vector_Length_Squared(result),  0.0001);
    }

    {
        // Rotate (0,y,0) 90 degrees around x-axis == (0,0,y)
        const double to_rotate[3]       = {0.0, 2.0, 0.0};
        const double rotation_axis[3]   = {1.0, 0.0, 0.0};  // z-axis
        const double rotation_angle     = 90.0 * deg_to_rad;
        double result[3];
        const double answer[3]          = {0.0, 0.0, to_rotate[1]};
        Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result);

        //std::cout << "to_rotate = (" << to_rotate[0] << ", " << to_rotate[1] << ", " << to_rotate[2] << ")\n";
        //std::cout << "answer = (" << answer[0] << ", " << answer[1] << ", " << answer[2] << ")\n";
        //std::cout << "result = (" << result[0] << ", " << result[1] << ", " << result[2] << ")\n";

        BOOST_CHECK_SMALL(result[0] - answer[0], 1.0e-7);
        BOOST_CHECK_SMALL(result[1] - answer[1], 1.0e-7);
        BOOST_CHECK_SMALL(result[2] - answer[2], 1.0e-7);

        BOOST_CHECK_CLOSE(2.0,                              Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length(to_rotate),         Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length_Squared(to_rotate), Vector_Length_Squared(result),  0.0001);
    }

    {
        // Rotate (0,y,0) 45 degrees around x-axis == sqrt(2)/2 * (0,y,y)
        const double to_rotate[3]       = {0.0, 2.0, 0.0};
        const double rotation_axis[3]   = {1.0, 0.0, 0.0};  // z-axis
        const double rotation_angle     = 45.0 * deg_to_rad;
        double result[3];
        const double answer[3]          = {0.0, 0.5*sqrt_2*to_rotate[1], 0.5*sqrt_2*to_rotate[1]};
        Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result);

        //std::cout << "to_rotate = (" << to_rotate[0] << ", " << to_rotate[1] << ", " << to_rotate[2] << ")\n";
        //std::cout << "answer = (" << answer[0] << ", " << answer[1] << ", " << answer[2] << ")\n";
        //std::cout << "result = (" << result[0] << ", " << result[1] << ", " << result[2] << ")\n";

        BOOST_CHECK_SMALL(result[0] - answer[0], 1.0e-7);
        BOOST_CHECK_SMALL(result[1] - answer[1], 1.0e-7);
        BOOST_CHECK_SMALL(result[2] - answer[2], 1.0e-7);

        BOOST_CHECK_CLOSE(2.0,                              Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length(to_rotate),         Vector_Length(result),          0.0001);
        BOOST_CHECK_CLOSE(Vector_Length_Squared(to_rotate), Vector_Length_Squared(result),  0.0001);
    }

    {
        // Throw an error if vector to rotate has length 0, only in debug build.
        const double to_rotate[3]       = {0.0, 0.0, 0.0};
        const double rotation_axis[3]   = {1.0, 0.0, 0.0};  // z-axis
        const double rotation_angle     = 45.0 * deg_to_rad;
        double result[3];

#ifdef YDEBUG
        BOOST_CHECK_THROW(Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result), std::runtime_error);
#endif // #ifdef YDEBUG
    }

    {
        // Throw an error if rotation axis is not unitary, only in debug build.
        const double to_rotate[3]       = {0.0, 2.0, 0.0};
        const double rotation_axis[3]   = {1.0, 2.0, 3.0};  // z-axis
        const double rotation_angle     = 45.0 * deg_to_rad;
        double result[3];

#ifdef YDEBUG
        BOOST_CHECK_THROW(Rodrigues_Rotation(to_rotate, rotation_axis, rotation_angle, result), std::runtime_error);
#endif // #ifdef YDEBUG
    }
}

BOOST_AUTO_TEST_SUITE_END()

