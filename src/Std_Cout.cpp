
#include <iomanip>

#include "Std_Cout.hpp"
// #include "Assert.hpp"


File_And_Screen_Stream std_cout;


// **************************************************************
File_And_Screen_Stream::~File_And_Screen_Stream(void)
/**
 * Destructor. Close the file if open.
 */
{
    Flush();

    if (filestream.is_open())
        filestream.close();
}

// **************************************************************
File_And_Screen_Stream & File_And_Screen_Stream::operator<<(std::ostream& (*pfun)(std::ostream&))
/**
 * This allow sending "std::endl" to the stream.
 */
{
    pfun(filestream);
    pfun(std::cout);
    return *this;
}

// **************************************************************
void File_And_Screen_Stream::open(std::string filename, const bool append)
/**
 * Open file.
 */
{
    std::cout << "Opening file " << filename << "...\n";
    if (append)
        filestream.open(filename.c_str(), std::ios_base::app);
    else
        filestream.open(filename.c_str(), std::ios_base::out);
    assert(filestream.is_open());
}

// **************************************************************
std::streamsize File_And_Screen_Stream::precision(const std::streamsize p)
/**
 * Set precision of stream.
 */
{
    filestream.precision(p);
    return std::cout.precision(p);
}

// **************************************************************
void File_And_Screen_Stream::Flush()
/**
 * Flush both file and screen output.
 */
{
    filestream << std::flush;
    std::cout  << std::flush;
}

// **************************************************************
void File_And_Screen_Stream::Format(const int width,
                                    const int nb_after_dot,
                                    const char type,
                                    const char justify,
                                    const char fill)
/**
 * Set the output format. Similar to printf() format.
 * @param width         Output's width (in characters)
 * @param nb_after_dot  Number of digits after decimal
 * @param type          Character describing the type
 * @param justify       Left (l) or right (r) justification.
 *                      Default to right.
 * @param fill          Filling character (default to space)
 *
 * Examples:
 *      Format(3, 0, 'd')       equivalent to printf("%3d", [...])
 *      Format(0, 4, 'f')       equivalent to printf("%.4f", [...])
 *      Format(11, 4, 'f')      equivalent to printf("%11.4f", [...])
 *      Format(12, 6, 'g', 'l') equivalent to printf("%-12.6g", [...])
 *      Format(12, 6, 'g', '-') equivalent to printf("%-12.6g", [...])
 *      Format(12, 6, 'g', '+') equivalent to printf("%+12.6g", [...])
 */
{
    if (width != 0)
    {
        std::cout  << std::setw(width);
        filestream << std::setw(width);
    }

    if (type != 'd')
    {
        // If not an integer, set the precision.
        std::cout  << std::setprecision(nb_after_dot);
        filestream << std::setprecision(nb_after_dot);

        if (type == 'f')
        {
            std::cout  << std::fixed;
            filestream << std::fixed;
        }
        else if (type == 'e')
        {
            std::cout  << std::scientific;
            filestream << std::scientific;
        }
    }
    if (justify == 'r')
    {
        std::cout  << std::right;
        filestream << std::right;
    } else if (justify == '+')
    {
        std::cout  << std::showpos;
        filestream << std::showpos;
    } else if (justify == 'l' || justify == '-')
    {
        std::cout  << std::left;
        filestream << std::left;
    }
    std::cout  << std::setfill(fill);
    filestream << std::setfill(fill);
}

// ********** End of file ***************************************