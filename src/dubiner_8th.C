#include "dubiner.h"

void Dubiner::dubiner_8th(const mpfr_class & zeta0,
                          const mpfr_class & zeta1,
                          const mpfr_class & zeta2,
                          std::vector<mpfr_class> & vals)
{
  /* (0,8) */ vals.push_back(9*pow(zeta2, 8) + 144*pow(zeta2, 7)*(2*zeta2 - 2) + 588*pow(zeta2, 6)*pow(2*zeta2 - 2, 2) + 882*pow(zeta2, 5)*pow(2*zeta2 - 2, 3) + (2205.0/4.0)*pow(zeta2, 4)*pow(2*zeta2 - 2, 4) + 147*pow(zeta2, 3)*pow(2*zeta2 - 2, 5) + (63.0/4.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 6) + (9.0/16.0)*zeta2*pow(2*zeta2 - 2, 7) + (1.0/256.0)*pow(2*zeta2 - 2, 8));
  /* (1,7) */ vals.push_back((-zeta0 + zeta1)*(120*pow(zeta2, 7) + 735*pow(zeta2, 6)*(2*zeta2 - 2) + 1323*pow(zeta2, 5)*pow(2*zeta2 - 2, 2) + (3675.0/4.0)*pow(zeta2, 4)*pow(2*zeta2 - 2, 3) + (525.0/2.0)*pow(zeta2, 3)*pow(2*zeta2 - 2, 4) + (945.0/32.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 5) + (35.0/32.0)*zeta2*pow(2*zeta2 - 2, 6) + (1.0/128.0)*pow(2*zeta2 - 2, 7)));
  /* (2,6) */ vals.push_back((pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2))*(462*pow(zeta2, 6) + 1386*pow(zeta2, 5)*(2*zeta2 - 2) + (2475.0/2.0)*pow(zeta2, 4)*pow(2*zeta2 - 2, 2) + (825.0/2.0)*pow(zeta2, 3)*pow(2*zeta2 - 2, 3) + (825.0/16.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 4) + (33.0/16.0)*zeta2*pow(2*zeta2 - 2, 5) + (1.0/64.0)*pow(2*zeta2 - 2, 6)));
  /* (3,5) */ vals.push_back((-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3))*(792*pow(zeta2, 5) + (2475.0/2.0)*pow(zeta2, 4)*(2*zeta2 - 2) + 550*pow(zeta2, 3)*pow(2*zeta2 - 2, 2) + (165.0/2.0)*pow(zeta2, 2)*pow(2*zeta2 - 2, 3) + (15.0/4.0)*zeta2*pow(2*zeta2 - 2, 4) + (1.0/32.0)*pow(2*zeta2 - 2, 5)));
  /* (4,4) */ vals.push_back((pow(zeta0, 4) - 16*pow(zeta0, 3)*zeta1 + 36*pow(zeta0, 2)*pow(zeta1, 2) - 16*zeta0*pow(zeta1, 3) + pow(zeta1, 4))*(715*pow(zeta2, 4) + 572*pow(zeta2, 3)*(2*zeta2 - 2) + 117*pow(zeta2, 2)*pow(2*zeta2 - 2, 2) + (13.0/2.0)*zeta2*pow(2*zeta2 - 2, 3) + (1.0/16.0)*pow(2*zeta2 - 2, 4)));
  /* (5,3) */ vals.push_back((364*pow(zeta2, 3) + (273.0/2.0)*pow(zeta2, 2)*(2*zeta2 - 2) + (21.0/2.0)*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3))*(-pow(zeta0, 5) + 25*pow(zeta0, 4)*zeta1 - 100*pow(zeta0, 3)*pow(zeta1, 2) + 100*pow(zeta0, 2)*pow(zeta1, 3) - 25*zeta0*pow(zeta1, 4) + pow(zeta1, 5)));
  /* (6,2) */ vals.push_back((105*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (15.0/2.0)*pow(2*zeta2 - 1, 2) - 15.0/2.0)*(pow(zeta0, 6) - 36*pow(zeta0, 5)*zeta1 + 225*pow(zeta0, 4)*pow(zeta1, 2) - 400*pow(zeta0, 3)*pow(zeta1, 3) + 225*pow(zeta0, 2)*pow(zeta1, 4) - 36*zeta0*pow(zeta1, 5) + pow(zeta1, 6)));
  /* (7,1) */ vals.push_back((17*zeta2 - 1)*(-pow(zeta0, 7) + 49*pow(zeta0, 6)*zeta1 - 441*pow(zeta0, 5)*pow(zeta1, 2) + 1225*pow(zeta0, 4)*pow(zeta1, 3) - 1225*pow(zeta0, 3)*pow(zeta1, 4) + 441*pow(zeta0, 2)*pow(zeta1, 5) - 49*zeta0*pow(zeta1, 6) + pow(zeta1, 7)));
  /* (8,0) */ vals.push_back(pow(zeta0, 8) - 64*pow(zeta0, 7)*zeta1 + 784*pow(zeta0, 6)*pow(zeta1, 2) - 3136*pow(zeta0, 5)*pow(zeta1, 3) + 4900*pow(zeta0, 4)*pow(zeta1, 4) - 3136*pow(zeta0, 3)*pow(zeta1, 5) + 784*pow(zeta0, 2)*pow(zeta1, 6) - 64*zeta0*pow(zeta1, 7) + pow(zeta1, 8));
}

// Local Variables:
// truncate-lines: t
// End:
