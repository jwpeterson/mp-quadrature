#include "dubiner.h"

void Dubiner::dubiner_5th(const mpfr_class & zeta0,
                          const mpfr_class & zeta1,
                          const mpfr_class & zeta2,
                          std::vector<mpfr_class> & vals)
{
  /* (0,5) */ vals.push_back(6*pow(zeta2, 5) + (75.0/2.0)*pow(zeta2, 4)*(2*zeta2 - 2) + 50*pow(zeta2, 3)*pow(2*zeta2 - 2, 2) + (75.0/4.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 3) + (15.0/8.0)*zeta2*pow(2*zeta2 - 2, 4) + (1.0/32.0)*pow(2*zeta2 - 2, 5));
  /* (1,4) */ vals.push_back((-zeta0 + zeta1)*(35*pow(zeta2, 4) + 70*pow(zeta2, 3)*(2*zeta2 - 2) + (63.0/2.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 2) + (7.0/2.0)*zeta2*pow(2*zeta2 - 2, 3) + (1.0/16.0)*pow(2*zeta2 - 2, 4)));
  /* (2,3) */ vals.push_back((pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2))*(56*pow(zeta2, 3) + 42*pow(zeta2, 2)*(2*zeta2 - 2) + 6*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3)));
  /* (3,2) */ vals.push_back((-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3))*(36*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (9.0/2.0)*pow(2*zeta2 - 1, 2) - 9.0/2.0));
  /* (4,1) */ vals.push_back((11*zeta2 - 1)*(pow(zeta0, 4) - 16*pow(zeta0, 3)*zeta1 + 36*pow(zeta0, 2)*pow(zeta1, 2) - 16*zeta0*pow(zeta1, 3) + pow(zeta1, 4)));
  /* (5,0) */ vals.push_back(-pow(zeta0, 5) + 25*pow(zeta0, 4)*zeta1 - 100*pow(zeta0, 3)*pow(zeta1, 2) + 100*pow(zeta0, 2)*pow(zeta1, 3) - 25*zeta0*pow(zeta1, 4) + pow(zeta1, 5));
}
