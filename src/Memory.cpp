
#include <stdint.h> // (u)int64_t

#include "Memory.hpp"

// **************************************************************
void Print_Double_in_Binary(double d)
/**
 * Prints binary representation of a double
 * http://www.exploringbinary.com/displaying-the-raw-fields-of-a-floating-point-number/
 */
{
    uint64_t *double_as_int = (uint64_t *) &d;
    const int bit_size = CHAR_BIT*sizeof(uint64_t);

    // Print bits by bits
    for (int b = 0 ; b <= bit_size-1 ; b++)
    {
        if (b == 1)
            std_cout << " ";    // Space after sign field
        if (b == 12)
            std_cout << " ";    // Space after exponent field

        // Get bit, but in reverse order. On Little Endian machines
        // (most of Intel and such), the byte with lower address
        // is the less significant. Since we want to print from
        // the most significant, we iterate from the end.
        if ((*double_as_int >> ((bit_size-1)-b)) & 1)
            std_cout << "1";
        else
            std_cout << "0";
    }
    //std_cout << "\n";
}


// ********** End of file ***************************************
