#ifndef INC_STRUCTURE_ELEMENT_HPP
#define INC_STRUCTURE_ELEMENT_HPP

#include "FloatType.hpp"

const int LASER_IS_PHOTONS  = 1;
const int LASER_IS_FIELD    = 2;
const int LASER_IS_BOTH     = 3;

class laser_props
{
   public:
   fdouble fwhm;
   fdouble ww;
   fdouble ww_au;
   fdouble la;
   fdouble Ep;
   fdouble intensity;
   fdouble pho_engy_au;
   fdouble nb_photons;
   bool cst_in_space;
   bool enable;
   fdouble k;
   fdouble field[3];
   int treated_as;

   bool treated_as_photons() { return (treated_as == LASER_IS_PHOTONS ? true : false); }
   bool treated_as_field()   { return (treated_as == LASER_IS_FIELD   ? true : false); }
   bool treated_as_both()    { return (treated_as == LASER_IS_BOTH    ? true : false); }
};

#endif // #ifndef INC_STRUCTURE_ELEMENT_HPP

// ********** End of file ***************************************
