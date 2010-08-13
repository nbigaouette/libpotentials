#ifndef INC_MEMORY_hpp
#define INC_MEMORY_hpp

#include <climits> // CHAR_BIT

#include "Std_Cout.hpp"

// **************************************************************
template <class Integer>
static inline void * calloc_and_check(Integer nb, size_t s, const std::string msg = "")
{
    void *p = NULL;
    p = calloc(nb, s);
    if (p == NULL)
    {
        double nb_s = double(nb) * double(s);
        std_cout << "Allocation of " << nb << " x " << s << " bytes = " << nb_s << " bytes failed" << std::endl;
        std_cout << "(" << double(nb_s) / (1024.0) << " KiB, "
                         << double(nb_s) / (1024.0*1024.0) << " MiB, "
                         << double(nb_s) / (1024.0*1024.0*1024.0) << " GiB)" << std::endl;
        std_cout << "p = " << p << "\n";
        if (msg != "")
        {
            std_cout << "Comment: " << msg << std::endl;
        }
        std_cout << "Aborting\n" << std::flush;
        abort();
    }
    return p;
}

// **************************************************************
template <class Integer>
static inline void * malloc_and_check(Integer nb, size_t s, const std::string msg = "")
{
    void *p = NULL;
    p = malloc(nb * s);
    if (p == NULL)
    {
        const long unsigned int slu  = (long unsigned int) s;
        const long unsigned int slun = (long unsigned int) s * nb;
        std_cout << "Allocation of " << nb << " x " << slu << " bytes = " << slun << " bytes failed" << std::endl;
        std_cout << "(" << double(slun) / (1024.0) << " KiB, "
                         << double(slun) / (1024.0*1024.0) << " MiB, "
                         << double(slun) / (1024.0*1024.0*1024.0) << " GiB)" << std::endl;
        std_cout << "p = " << p << "\n";
        if (msg != "")
        {
            std_cout << "Comment: " << msg << std::endl;
        }
        std_cout << "Aborting\n" << std::flush;
        abort();
    }
    return p;
}

// **************************************************************
template <class Pointer>
void free_me(Pointer &p)
{
    if (p != NULL)
    {
        free(p);
    }
    p = NULL;
}

void Print_Double_in_Binary(double d);

// **************************************************************
template <class Integer>
void Print_Integer_in_Binary(Integer n)
/**
 * Prints binary representation of an integer if any size.
 * Inspired by http://www.exploringbinary.com/displaying-the-raw-fields-of-a-floating-point-number/
 * WARNING: In C/C++, logical right shift of SIGNED integers is compiler dependant. GCC keeps the
 *          sign bit intact (instead of putting a 0).
 *          So ">>" is an arithmetic shift when the integer is signed. Unsigned are not
 *          affected (arithmetic and logical shifts are the same for unsigned integers).
 *          See http://en.wikipedia.org/wiki/Bitwise_operation#Arithmetic_shift
 */
{
                                        // Example 32 bits integers, converted from
                                        // http://www.binaryconvert.com/convert_unsigned_int.html
    const Integer zero  =  Integer(0);  // 00000000 00000000 00000000 00000000
    //const Integer ones  = ~zero;        // 11111111 11111111 11111111 11111111
    const Integer one   =  Integer(1);  // 00000000 00000000 00000000 00000001
    //const Integer two   =  Integer(2);  // 00000000 00000000 00000000 00000010
    //const Integer eigth =  Integer(8);  // 00000000 00000000 00000000 00001000
    const Integer nb_bits_per_byte    = CHAR_BIT; // Normaly, it is 8, but could be different.
    const Integer nb_bits_per_Integer = sizeof(n)*nb_bits_per_byte;

    // Starting from the LSB being index "0", the MSB is at index "msb_position"
    const Integer msb_position  = nb_bits_per_Integer - one;
    const Integer msb           = one << msb_position;
    const Integer or_msb        = ~msb;

    // Note that right shifting a signed integer migth keep the sign bit intact
    // (instead of setting it to 0) because C/C++ is implementation dependant
    // regarding right shift applied to negative signed integers. GCC will do
    // an "arithmetic right shift", meaning dividing the integer by 2. This will
    // keep the number negative (if it was). Because of this, the mask can get
    // screwed. If the Integer type is signed, first right shifting of the
    // mask of one (having an initial value of "msb" == 10000... and thus a
    // negative value) will keep the sign bit (leading to mask == 11000...) but
    // what we want is just to move the mask's bit, not keep the integer
    // reprentation "valid" (we want mask == 01000...). To fix that, after
    // right shifting the mask by one, we "AND" it (using "&") with "or_msb"
    // (or_msb == 01111...) to make sure we forget the sign bit.
    for (Integer mask = msb ; mask != zero ; mask = ((mask >> one) & or_msb ))
    {
        // If "n"'s bit at position of the mask is 0, print 0, else print 1.
        if ((mask & n) == zero) std_cout << "0";
        else                    std_cout << "1";
    }
}


#endif // INC_MEMORY_hpp

// ********** End of file ***************************************
