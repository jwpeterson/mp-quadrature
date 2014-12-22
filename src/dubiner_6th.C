#include "dubiner.h"

void Dubiner::dubiner_6th(const mpfr_class & zeta0,
                          const mpfr_class & zeta1,
                          const mpfr_class & zeta2,
                          std::vector<mpfr_class> & vals)
{
  /* (0,6) */ vals.push_back(7*pow(zeta2, 6) + 63*pow(zeta2, 5)*(2*zeta2 - 2) + (525.0/4.0)*pow(zeta2, 4)*pow(2*zeta2 - 2, 2) + (175.0/2.0)*pow(zeta2, 3)*pow(2*zeta2 - 2, 3) + (315.0/16.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 4) + (21.0/16.0)*zeta2*pow(2*zeta2 - 2, 5) + (1.0/64.0)*pow(2*zeta2 - 2, 6));
  /* (1,5) */ vals.push_back((-zeta0 + zeta1)*(56*pow(zeta2, 5) + 175*pow(zeta2, 4)*(2*zeta2 - 2) + 140*pow(zeta2, 3)*pow(2*zeta2 - 2, 2) + 35*pow(zeta2, 2)*pow(2*zeta2 - 2, 3) + (5.0/2.0)*zeta2*pow(2*zeta2 - 2, 4) + (1.0/32.0)*pow(2*zeta2 - 2, 5)));
  /* (2,4) */ vals.push_back((pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2))*(126*pow(zeta2, 4) + 168*pow(zeta2, 3)*(2*zeta2 - 2) + 54*pow(zeta2, 2)*pow(2*zeta2 - 2, 2) + (9.0/2.0)*zeta2*pow(2*zeta2 - 2, 3) + (1.0/16.0)*pow(2*zeta2 - 2, 4)));
  /* (3,3) */ vals.push_back((-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3))*(120*pow(zeta2, 3) + (135.0/2.0)*pow(zeta2, 2)*(2*zeta2 - 2) + (15.0/2.0)*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3)));
  /* (4,2) */ vals.push_back((55*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (11.0/2.0)*pow(2*zeta2 - 1, 2) - 11.0/2.0)*(pow(zeta0, 4) - 16*pow(zeta0, 3)*zeta1 + 36*pow(zeta0, 2)*pow(zeta1, 2) - 16*zeta0*pow(zeta1, 3) + pow(zeta1, 4)));
  /* (5,1) */ vals.push_back((13*zeta2 - 1)*(-pow(zeta0, 5) + 25*pow(zeta0, 4)*zeta1 - 100*pow(zeta0, 3)*pow(zeta1, 2) + 100*pow(zeta0, 2)*pow(zeta1, 3) - 25*zeta0*pow(zeta1, 4) + pow(zeta1, 5)));
  /* (6,0) */ vals.push_back(pow(zeta0, 6) - 36*pow(zeta0, 5)*zeta1 + 225*pow(zeta0, 4)*pow(zeta1, 2) - 400*pow(zeta0, 3)*pow(zeta1, 3) + 225*pow(zeta0, 2)*pow(zeta1, 4) - 36*zeta0*pow(zeta1, 5) + pow(zeta1, 6));
}

// Local Variables:
// truncate-lines: t
// End:
