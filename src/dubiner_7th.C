#include "dubiner.h"

void Dubiner::dubiner_7th(const mpfr_class & zeta0,
                          const mpfr_class & zeta1,
                          const mpfr_class & zeta2,
                          std::vector<mpfr_class> & vals)
{
  /* (0,7) */ vals.push_back(8*pow(zeta2, 7) + 98*pow(zeta2, 6)*(2*zeta2 - 2) + 294*pow(zeta2, 5)*pow(2*zeta2 - 2, 2) + (1225.0/4.0)*pow(zeta2, 4)*pow(2*zeta2 - 2, 3) + (245.0/2.0)*pow(zeta2, 3)*pow(2*zeta2 - 2, 4) + (147.0/8.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 5) + (7.0/8.0)*zeta2*pow(2*zeta2 - 2, 6) + (1.0/128.0)*pow(2*zeta2 - 2, 7));
  /* (1,6) */ vals.push_back((-zeta0 + zeta1)*(84*pow(zeta2, 6) + 378*pow(zeta2, 5)*(2*zeta2 - 2) + (945.0/2.0)*pow(zeta2, 4)*pow(2*zeta2 - 2, 2) + 210*pow(zeta2, 3)*pow(2*zeta2 - 2, 3) + (135.0/4.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 4) + (27.0/16.0)*zeta2*pow(2*zeta2 - 2, 5) + (1.0/64.0)*pow(2*zeta2 - 2, 6)));
  /* (2,5) */ vals.push_back((pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2))*(252*pow(zeta2, 5) + 525*pow(zeta2, 4)*(2*zeta2 - 2) + 300*pow(zeta2, 3)*pow(2*zeta2 - 2, 2) + (225.0/4.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 3) + (25.0/8.0)*zeta2*pow(2*zeta2 - 2, 4) + (1.0/32.0)*pow(2*zeta2 - 2, 5)));
  /* (3,4) */ vals.push_back((-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3))*(330*pow(zeta2, 4) + 330*pow(zeta2, 3)*(2*zeta2 - 2) + (165.0/2.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 2) + (11.0/2.0)*zeta2*pow(2*zeta2 - 2, 3) + (1.0/16.0)*pow(2*zeta2 - 2, 4)));
  /* (4,3) */ vals.push_back((220*pow(zeta2, 3) + 99*pow(zeta2, 2)*(2*zeta2 - 2) + 9*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3))*(pow(zeta0, 4) - 16*pow(zeta0, 3)*zeta1 + 36*pow(zeta0, 2)*pow(zeta1, 2) - 16*zeta0*pow(zeta1, 3) + pow(zeta1, 4)));
  /* (5,2) */ vals.push_back((78*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (13.0/2.0)*pow(2*zeta2 - 1, 2) - 13.0/2.0)*(-pow(zeta0, 5) + 25*pow(zeta0, 4)*zeta1 - 100*pow(zeta0, 3)*pow(zeta1, 2) + 100*pow(zeta0, 2)*pow(zeta1, 3) - 25*zeta0*pow(zeta1, 4) + pow(zeta1, 5)));
  /* (6,1) */ vals.push_back((15*zeta2 - 1)*(pow(zeta0, 6) - 36*pow(zeta0, 5)*zeta1 + 225*pow(zeta0, 4)*pow(zeta1, 2) - 400*pow(zeta0, 3)*pow(zeta1, 3) + 225*pow(zeta0, 2)*pow(zeta1, 4) - 36*zeta0*pow(zeta1, 5) + pow(zeta1, 6)));
  /* (7,0) */ vals.push_back(-pow(zeta0, 7) + 49*pow(zeta0, 6)*zeta1 - 441*pow(zeta0, 5)*pow(zeta1, 2) + 1225*pow(zeta0, 4)*pow(zeta1, 3) - 1225*pow(zeta0, 3)*pow(zeta1, 4) + 441*pow(zeta0, 2)*pow(zeta1, 5) - 49*zeta0*pow(zeta1, 6) + pow(zeta1, 7));
}

// Local Variables:
// truncate-lines: t
// End:
