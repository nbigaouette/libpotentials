#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE LibPotentialsTesting
#include <boost/test/unit_test.hpp>

/**
 * To add a new test, create a new .cpp file.
 * You can inspire yourself from "unit_testing/example.cpp".
 * Note that "BOOST_TEST_MODULE" can only be defined once (in
 * this case here) and should not be re-defined in other files.
 *
 * This is based on Boost's Test Library.
 * http://www.boost.org/doc/libs/1_45_0/libs/test/doc/html/index.html
 * You'll need to install boost.
 */

// /*

// See http://bitten.edgewall.org/wiki/BoostTest
#include <iostream>
#include <fstream>
// #include <cassert>


// Redirect code's output (std::cout) to log file
struct ReportRedirector
{
    std::streambuf *orig;
    std::ofstream out;

    ReportRedirector() : out("output/unit_test.log")
    {
        assert(out.is_open());
        orig = std::cout.rdbuf(out.rdbuf());
    }

    ~ReportRedirector()
    {
        std::cout.rdbuf(orig);
    }
};

static ReportRedirector foo;

// */
