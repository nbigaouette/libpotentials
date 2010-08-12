#ifndef INC_STRUCTURE_ELEMENT_HPP
#define INC_STRUCTURE_ELEMENT_HPP

const int LASER_IS_PHOTONS  = 1;
const int LASER_IS_FIELD    = 2;
const int LASER_IS_BOTH     = 3;

class laser_props
{
   public:
   double fwhm;
   double ww;
   double ww_au;
   double la;
   double Ep;
   double intensity;
   double pho_engy_au;
   double nb_photons;
   bool cst_in_space;
   bool enable;
   double k;
   double field[3];
   int treated_as;

   bool treated_as_photons() { return (treated_as == LASER_IS_PHOTONS ? true : false); }
   bool treated_as_field()   { return (treated_as == LASER_IS_FIELD   ? true : false); }
   bool treated_as_both()    { return (treated_as == LASER_IS_BOTH    ? true : false); }
};

#endif // #ifndef INC_STRUCTURE_ELEMENT_HPP

// ********** End of file ***************************************
