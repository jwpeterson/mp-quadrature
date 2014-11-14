#include "dubiner.h"

void Dubiner::p(unsigned d,
                const mpfr_class & xi,
                const mpfr_class & eta,
                std::vector<mpfr_class> & vals)
{
  // Make sure there are no old values in the vector.
  vals.clear();

  // Compute the barycenter coordinate values
  mpfr_class
    zeta0 = xi,
    zeta1 = eta,
    zeta2 = 1. - xi - eta;

  /* (0,0) */ vals.push_back(1);

  if (d >= 1)
    {
      /* (0,1) */ vals.push_back(3*zeta2 - 1);
      /* (1,0) */ vals.push_back(-zeta0 + zeta1);
    }
  if (d >= 2)
    {
      // Note: gmpfrxx does not like it when you have an "L" suffix on your constants!
      // We need to be careful here if there are fractions that do not have a finite
      // decimal representation -- they will not be full precision unless you explicitly
      // cast stuff to the mpfr_class type...
      /* (0,2) */ vals.push_back(3*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (3.0/2.0)*pow(2*zeta2 - 1, 2) - 3.0/2.0);
      /* (1,1) */ vals.push_back((-zeta0 + zeta1)*(5*zeta2 - 1));
      /* (2,0) */ vals.push_back(pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2));
    }
  if (d >= 3)
    {
      /* (0,3) */ vals.push_back(4*pow(zeta2, 3) + 9*pow(zeta2, 2)*(2*zeta2 - 2) + 3*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3));
      /* (1,2) */ vals.push_back((-zeta0 + zeta1)*(10*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (5.0/2.0)*pow(2*zeta2 - 1, 2) - 5.0/2.0));
      /* (2,1) */ vals.push_back((7*zeta2 - 1)*(pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2)));
      /* (3,0) */ vals.push_back(-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3));
    }
  if (d >= 4)
    {
      /* (0,4) */ vals.push_back(5*pow(zeta2, 4) + 20*pow(zeta2, 3)*(2*zeta2 - 2) + 15*pow(zeta2, 2)*pow(2*zeta2 - 2, 2) + (5.0/2.0)*zeta2*pow(2*zeta2 - 2, 3) + (1.0/16.0)*pow(2*zeta2 - 2, 4));
      /* (1,3) */ vals.push_back((-zeta0 + zeta1)*(20*pow(zeta2, 3) + (45.0/2.0)*pow(zeta2, 2)*(2*zeta2 - 2) + (9.0/2.0)*zeta2*pow(2*zeta2 - 2, 2) + (1.0/8.0)*pow(2*zeta2 - 2, 3)));
      /* (2,2) */ vals.push_back((pow(zeta0, 2) - 4*zeta0*zeta1 + pow(zeta1, 2))*(21*pow(zeta2, 2) + (1.0/4.0)*pow(2*zeta2 - 2, 2) + (7.0/2.0)*pow(2*zeta2 - 1, 2) - 7.0/2.0));
      /* (3,1) */ vals.push_back((9*zeta2 - 1)*(-pow(zeta0, 3) + 9*pow(zeta0, 2)*zeta1 - 9*zeta0*pow(zeta1, 2) + pow(zeta1, 3)));
      /* (4,0) */ vals.push_back(pow(zeta0, 4) - 16*pow(zeta0, 3)*zeta1 + 36*pow(zeta0, 2)*pow(zeta1, 2) - 16*zeta0*pow(zeta1, 3) + pow(zeta1, 4));
    }
  if (d >= 5)
    this->dubiner_5th(zeta0, zeta1, zeta2, vals);
  if (d >= 6)
    this->dubiner_6th(zeta0, zeta1, zeta2, vals);
  if (d >= 7)
    this->dubiner_7th(zeta0, zeta1, zeta2, vals);
  if (d >= 8)
    this->dubiner_8th(zeta0, zeta1, zeta2, vals);
  if (d >= 9)
    this->dubiner_9th(zeta0, zeta1, zeta2, vals);
  if (d >= 10)
    this->dubiner_10th(zeta0, zeta1, zeta2, vals);
  if (d >= 11)
    this->dubiner_11th(zeta0, zeta1, zeta2, vals);
  if (d >= 12)
    this->dubiner_12th(zeta0, zeta1, zeta2, vals);
  if (d >= 13)
    this->dubiner_13th(zeta0, zeta1, zeta2, vals);
  if (d >= 14)
    this->dubiner_14th(zeta0, zeta1, zeta2, vals);
  if (d >= 15)
    this->dubiner_15th(zeta0, zeta1, zeta2, vals);
  if (d >= 16)
    this->dubiner_16th(zeta0, zeta1, zeta2, vals);
  if (d >= 17)
    this->dubiner_17th(zeta0, zeta1, zeta2, vals);
  if (d >= 18)
    this->dubiner_18th(zeta0, zeta1, zeta2, vals);
  if (d >= 19)
    this->dubiner_19th(zeta0, zeta1, zeta2, vals);
  if (d >= 20)
    this->dubiner_20th(zeta0, zeta1, zeta2, vals);
  if (d >= 21)
    this->dubiner_21st(zeta0, zeta1, zeta2, vals);
  if (d >= 22)
    this->dubiner_22nd(zeta0, zeta1, zeta2, vals);
  if (d >= 23)
    this->dubiner_23rd(zeta0, zeta1, zeta2, vals);
  if (d >= 24)
    this->dubiner_24th(zeta0, zeta1, zeta2, vals);
  if (d >= 25)
    this->dubiner_25th(zeta0, zeta1, zeta2, vals);
  if (d >= 26)
    this->dubiner_26th(zeta0, zeta1, zeta2, vals);
  if (d >= 27)
    this->dubiner_27th(zeta0, zeta1, zeta2, vals);
  if (d >= 28)
    this->dubiner_28th(zeta0, zeta1, zeta2, vals);
  if (d >= 29)
    this->dubiner_29th(zeta0, zeta1, zeta2, vals);
  if (d >= 30)
    this->dubiner_30th(zeta0, zeta1, zeta2, vals);
}


void Dubiner::orthogonality_coeffs(unsigned d,
                                   std::vector<mpq_class> & coeffs)
{
  coeffs.clear();

  coeffs.push_back(mpq_class("1/2"));

  if (d >= 1)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/4"));
      /* 1 */ coeffs.push_back(mpq_class("1/12"));
    }
  if (d >= 2)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/6"));
      /* 1 */ coeffs.push_back(mpq_class("1/18"));
      /* 2 */ coeffs.push_back(mpq_class("1/30"));
    }
  if (d >= 3)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/8"));
      /* 1 */ coeffs.push_back(mpq_class("1/24"));
      /* 2 */ coeffs.push_back(mpq_class("1/40"));
      /* 3 */ coeffs.push_back(mpq_class("1/56"));
    }
  if (d >= 4)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/10"));
      /* 1 */ coeffs.push_back(mpq_class("1/30"));
      /* 2 */ coeffs.push_back(mpq_class("1/50"));
      /* 3 */ coeffs.push_back(mpq_class("1/70"));
      /* 4 */ coeffs.push_back(mpq_class("1/90"));
    }
  if (d >= 5)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/12"));
      /* 1 */ coeffs.push_back(mpq_class("1/36"));
      /* 2 */ coeffs.push_back(mpq_class("1/60"));
      /* 3 */ coeffs.push_back(mpq_class("1/84"));
      /* 4 */ coeffs.push_back(mpq_class("1/108"));
      /* 5 */ coeffs.push_back(mpq_class("1/132"));
    }
  if (d >= 6)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/14"));
      /* 1 */ coeffs.push_back(mpq_class("1/42"));
      /* 2 */ coeffs.push_back(mpq_class("1/70"));
      /* 3 */ coeffs.push_back(mpq_class("1/98"));
      /* 4 */ coeffs.push_back(mpq_class("1/126"));
      /* 5 */ coeffs.push_back(mpq_class("1/154"));
      /* 6 */ coeffs.push_back(mpq_class("1/182"));
    }
  if (d >= 7)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/16"));
      /* 1 */ coeffs.push_back(mpq_class("1/48"));
      /* 2 */ coeffs.push_back(mpq_class("1/80"));
      /* 3 */ coeffs.push_back(mpq_class("1/112"));
      /* 4 */ coeffs.push_back(mpq_class("1/144"));
      /* 5 */ coeffs.push_back(mpq_class("1/176"));
      /* 6 */ coeffs.push_back(mpq_class("1/208"));
      /* 7 */ coeffs.push_back(mpq_class("1/240"));
    }
  if (d >= 8)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/18"));
      /* 1 */ coeffs.push_back(mpq_class("1/54"));
      /* 2 */ coeffs.push_back(mpq_class("1/90"));
      /* 3 */ coeffs.push_back(mpq_class("1/126"));
      /* 4 */ coeffs.push_back(mpq_class("1/162"));
      /* 5 */ coeffs.push_back(mpq_class("1/198"));
      /* 6 */ coeffs.push_back(mpq_class("1/234"));
      /* 7 */ coeffs.push_back(mpq_class("1/270"));
      /* 8 */ coeffs.push_back(mpq_class("1/306"));
    }
  if (d >= 9)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/20"));
      /* 1 */ coeffs.push_back(mpq_class("1/60"));
      /* 2 */ coeffs.push_back(mpq_class("1/100"));
      /* 3 */ coeffs.push_back(mpq_class("1/140"));
      /* 4 */ coeffs.push_back(mpq_class("1/180"));
      /* 5 */ coeffs.push_back(mpq_class("1/220"));
      /* 6 */ coeffs.push_back(mpq_class("1/260"));
      /* 7 */ coeffs.push_back(mpq_class("1/300"));
      /* 8 */ coeffs.push_back(mpq_class("1/340"));
      /* 9 */ coeffs.push_back(mpq_class("1/380"));
    }
  if (d >= 10)
    {
      /* 0 */  coeffs.push_back(mpq_class("1/22"));
      /* 1 */  coeffs.push_back(mpq_class("1/66"));
      /* 2 */  coeffs.push_back(mpq_class("1/110"));
      /* 3 */  coeffs.push_back(mpq_class("1/154"));
      /* 4 */  coeffs.push_back(mpq_class("1/198"));
      /* 5 */  coeffs.push_back(mpq_class("1/242"));
      /* 6 */  coeffs.push_back(mpq_class("1/286"));
      /* 7 */  coeffs.push_back(mpq_class("1/330"));
      /* 8 */  coeffs.push_back(mpq_class("1/374"));
      /* 9 */  coeffs.push_back(mpq_class("1/418"));
      /* 10 */ coeffs.push_back(mpq_class("1/462"));
    }
  if (d >= 11)
    {
      /* 0 */  coeffs.push_back(mpq_class("1/24"));
      /* 1 */  coeffs.push_back(mpq_class("1/72"));
      /* 2 */  coeffs.push_back(mpq_class("1/120"));
      /* 3 */  coeffs.push_back(mpq_class("1/168"));
      /* 4 */  coeffs.push_back(mpq_class("1/216"));
      /* 5 */  coeffs.push_back(mpq_class("1/264"));
      /* 6 */  coeffs.push_back(mpq_class("1/312"));
      /* 7 */  coeffs.push_back(mpq_class("1/360"));
      /* 8 */  coeffs.push_back(mpq_class("1/408"));
      /* 9 */  coeffs.push_back(mpq_class("1/456"));
      /* 10 */ coeffs.push_back(mpq_class("1/504"));
      /* 11 */ coeffs.push_back(mpq_class("1/552"));
    }
  if (d >= 12)
    {
      /* 0 */  coeffs.push_back(mpq_class("1/26"));
      /* 1 */  coeffs.push_back(mpq_class("1/78"));
      /* 2 */  coeffs.push_back(mpq_class("1/130"));
      /* 3 */  coeffs.push_back(mpq_class("1/182"));
      /* 4 */  coeffs.push_back(mpq_class("1/234"));
      /* 5 */  coeffs.push_back(mpq_class("1/286"));
      /* 6 */  coeffs.push_back(mpq_class("1/338"));
      /* 7 */  coeffs.push_back(mpq_class("1/390"));
      /* 8 */  coeffs.push_back(mpq_class("1/442"));
      /* 9 */  coeffs.push_back(mpq_class("1/494"));
      /* 10 */ coeffs.push_back(mpq_class("1/546"));
      /* 11 */ coeffs.push_back(mpq_class("1/598"));
      /* 12 */ coeffs.push_back(mpq_class("1/650"));
    }
  if (d >= 13)
    {
      /*  0 */ coeffs.push_back(mpq_class("1/28"));
      /*  1 */ coeffs.push_back(mpq_class("1/84"));
      /*  2 */ coeffs.push_back(mpq_class("1/140"));
      /*  3 */ coeffs.push_back(mpq_class("1/196"));
      /*  4 */ coeffs.push_back(mpq_class("1/252"));
      /*  5 */ coeffs.push_back(mpq_class("1/308"));
      /*  6 */ coeffs.push_back(mpq_class("1/364"));
      /*  7 */ coeffs.push_back(mpq_class("1/420"));
      /*  8 */ coeffs.push_back(mpq_class("1/476"));
      /*  9 */ coeffs.push_back(mpq_class("1/532"));
      /* 10 */ coeffs.push_back(mpq_class("1/588"));
      /* 11 */ coeffs.push_back(mpq_class("1/644"));
      /* 12 */ coeffs.push_back(mpq_class("1/700"));
      /* 13 */ coeffs.push_back(mpq_class("1/756"));
    }
  if (d >= 14)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/30"));
      /* 1  */ coeffs.push_back(mpq_class("1/90"));
      /* 2  */ coeffs.push_back(mpq_class("1/150"));
      /* 3  */ coeffs.push_back(mpq_class("1/210"));
      /* 4  */ coeffs.push_back(mpq_class("1/270"));
      /* 5  */ coeffs.push_back(mpq_class("1/330"));
      /* 6  */ coeffs.push_back(mpq_class("1/390"));
      /* 7  */ coeffs.push_back(mpq_class("1/450"));
      /* 8  */ coeffs.push_back(mpq_class("1/510"));
      /* 9  */ coeffs.push_back(mpq_class("1/570"));
      /* 10 */ coeffs.push_back(mpq_class("1/630"));
      /* 11 */ coeffs.push_back(mpq_class("1/690"));
      /* 12 */ coeffs.push_back(mpq_class("1/750"));
      /* 13 */ coeffs.push_back(mpq_class("1/810"));
      /* 14 */ coeffs.push_back(mpq_class("1/870"));
    }
  if (d >= 15)
    {
      /* 0 */  coeffs.push_back(mpq_class("1/32"));
      /* 1 */  coeffs.push_back(mpq_class("1/96"));
      /* 2 */  coeffs.push_back(mpq_class("1/160"));
      /* 3 */  coeffs.push_back(mpq_class("1/224"));
      /* 4 */  coeffs.push_back(mpq_class("1/288"));
      /* 5 */  coeffs.push_back(mpq_class("1/352"));
      /* 6 */  coeffs.push_back(mpq_class("1/416"));
      /* 7 */  coeffs.push_back(mpq_class("1/480"));
      /* 8 */  coeffs.push_back(mpq_class("1/544"));
      /* 9 */  coeffs.push_back(mpq_class("1/608"));
      /* 10 */ coeffs.push_back(mpq_class("1/672"));
      /* 11 */ coeffs.push_back(mpq_class("1/736"));
      /* 12 */ coeffs.push_back(mpq_class("1/800"));
      /* 13 */ coeffs.push_back(mpq_class("1/864"));
      /* 14 */ coeffs.push_back(mpq_class("1/928"));
      /* 15 */ coeffs.push_back(mpq_class("1/992"));
    }
  if (d >= 16)
    {
      /* 0 */  coeffs.push_back(mpq_class("1/34"));
      /* 1 */  coeffs.push_back(mpq_class("1/102"));
      /* 2 */  coeffs.push_back(mpq_class("1/170"));
      /* 3 */  coeffs.push_back(mpq_class("1/238"));
      /* 4 */  coeffs.push_back(mpq_class("1/306"));
      /* 5 */  coeffs.push_back(mpq_class("1/374"));
      /* 6 */  coeffs.push_back(mpq_class("1/442"));
      /* 7 */  coeffs.push_back(mpq_class("1/510"));
      /* 8 */  coeffs.push_back(mpq_class("1/578"));
      /* 9 */  coeffs.push_back(mpq_class("1/646"));
      /* 10 */ coeffs.push_back(mpq_class("1/714"));
      /* 11 */ coeffs.push_back(mpq_class("1/782"));
      /* 12 */ coeffs.push_back(mpq_class("1/850"));
      /* 13 */ coeffs.push_back(mpq_class("1/918"));
      /* 14 */ coeffs.push_back(mpq_class("1/986"));
      /* 15 */ coeffs.push_back(mpq_class("1/1054"));
      /* 16 */ coeffs.push_back(mpq_class("1/1122"));
    }
  if (d >= 17)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/36"));
      /* 1  */ coeffs.push_back(mpq_class("1/108"));
      /* 2  */ coeffs.push_back(mpq_class("1/180"));
      /* 3  */ coeffs.push_back(mpq_class("1/252"));
      /* 4  */ coeffs.push_back(mpq_class("1/324"));
      /* 5  */ coeffs.push_back(mpq_class("1/396"));
      /* 6  */ coeffs.push_back(mpq_class("1/468"));
      /* 7  */ coeffs.push_back(mpq_class("1/540"));
      /* 8  */ coeffs.push_back(mpq_class("1/612"));
      /* 9  */ coeffs.push_back(mpq_class("1/684"));
      /* 10 */ coeffs.push_back(mpq_class("1/756"));
      /* 11 */ coeffs.push_back(mpq_class("1/828"));
      /* 12 */ coeffs.push_back(mpq_class("1/900"));
      /* 13 */ coeffs.push_back(mpq_class("1/972"));
      /* 14 */ coeffs.push_back(mpq_class("1/1044"));
      /* 15 */ coeffs.push_back(mpq_class("1/1116"));
      /* 16 */ coeffs.push_back(mpq_class("1/1188"));
      /* 17 */ coeffs.push_back(mpq_class("1/1260"));
    }
  if (d >= 18)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/38"));
      /* 1  */ coeffs.push_back(mpq_class("1/114"));
      /* 2  */ coeffs.push_back(mpq_class("1/190"));
      /* 3  */ coeffs.push_back(mpq_class("1/266"));
      /* 4  */ coeffs.push_back(mpq_class("1/342"));
      /* 5  */ coeffs.push_back(mpq_class("1/418"));
      /* 6  */ coeffs.push_back(mpq_class("1/494"));
      /* 7  */ coeffs.push_back(mpq_class("1/570"));
      /* 8  */ coeffs.push_back(mpq_class("1/646"));
      /* 9  */ coeffs.push_back(mpq_class("1/722"));
      /* 10 */ coeffs.push_back(mpq_class("1/798"));
      /* 11 */ coeffs.push_back(mpq_class("1/874"));
      /* 12 */ coeffs.push_back(mpq_class("1/950"));
      /* 13 */ coeffs.push_back(mpq_class("1/1026"));
      /* 14 */ coeffs.push_back(mpq_class("1/1102"));
      /* 15 */ coeffs.push_back(mpq_class("1/1178"));
      /* 16 */ coeffs.push_back(mpq_class("1/1254"));
      /* 17 */ coeffs.push_back(mpq_class("1/1330"));
      /* 18 */ coeffs.push_back(mpq_class("1/1406"));
    }
  if (d >= 19)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/40"));
      /* 1  */ coeffs.push_back(mpq_class("1/120"));
      /* 2  */ coeffs.push_back(mpq_class("1/200"));
      /* 3  */ coeffs.push_back(mpq_class("1/280"));
      /* 4  */ coeffs.push_back(mpq_class("1/360"));
      /* 5  */ coeffs.push_back(mpq_class("1/440"));
      /* 6  */ coeffs.push_back(mpq_class("1/520"));
      /* 7  */ coeffs.push_back(mpq_class("1/600"));
      /* 8  */ coeffs.push_back(mpq_class("1/680"));
      /* 9  */ coeffs.push_back(mpq_class("1/760"));
      /* 10 */ coeffs.push_back(mpq_class("1/840"));
      /* 11 */ coeffs.push_back(mpq_class("1/920"));
      /* 12 */ coeffs.push_back(mpq_class("1/1000"));
      /* 13 */ coeffs.push_back(mpq_class("1/1080"));
      /* 14 */ coeffs.push_back(mpq_class("1/1160"));
      /* 15 */ coeffs.push_back(mpq_class("1/1240"));
      /* 16 */ coeffs.push_back(mpq_class("1/1320"));
      /* 17 */ coeffs.push_back(mpq_class("1/1400"));
      /* 18 */ coeffs.push_back(mpq_class("1/1480"));
      /* 19 */ coeffs.push_back(mpq_class("1/1560"));
    }
  if (d >= 20)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/42"));
      /* 1  */ coeffs.push_back(mpq_class("1/126"));
      /* 2  */ coeffs.push_back(mpq_class("1/210"));
      /* 3  */ coeffs.push_back(mpq_class("1/294"));
      /* 4  */ coeffs.push_back(mpq_class("1/378"));
      /* 5  */ coeffs.push_back(mpq_class("1/462"));
      /* 6  */ coeffs.push_back(mpq_class("1/546"));
      /* 7  */ coeffs.push_back(mpq_class("1/630"));
      /* 8  */ coeffs.push_back(mpq_class("1/714"));
      /* 9  */ coeffs.push_back(mpq_class("1/798"));
      /* 10 */ coeffs.push_back(mpq_class("1/882"));
      /* 11 */ coeffs.push_back(mpq_class("1/966"));
      /* 12 */ coeffs.push_back(mpq_class("1/1050"));
      /* 13 */ coeffs.push_back(mpq_class("1/1134"));
      /* 14 */ coeffs.push_back(mpq_class("1/1218"));
      /* 15 */ coeffs.push_back(mpq_class("1/1302"));
      /* 16 */ coeffs.push_back(mpq_class("1/1386"));
      /* 17 */ coeffs.push_back(mpq_class("1/1470"));
      /* 18 */ coeffs.push_back(mpq_class("1/1554"));
      /* 19 */ coeffs.push_back(mpq_class("1/1638"));
      /* 20 */ coeffs.push_back(mpq_class("1/1722"));
    }
  if (d >= 21)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/44"));
      /* 1  */ coeffs.push_back(mpq_class("1/132"));
      /* 2  */ coeffs.push_back(mpq_class("1/220"));
      /* 3  */ coeffs.push_back(mpq_class("1/308"));
      /* 4  */ coeffs.push_back(mpq_class("1/396"));
      /* 5  */ coeffs.push_back(mpq_class("1/484"));
      /* 6  */ coeffs.push_back(mpq_class("1/572"));
      /* 7  */ coeffs.push_back(mpq_class("1/660"));
      /* 8  */ coeffs.push_back(mpq_class("1/748"));
      /* 9  */ coeffs.push_back(mpq_class("1/836"));
      /* 10 */ coeffs.push_back(mpq_class("1/924"));
      /* 11 */ coeffs.push_back(mpq_class("1/1012"));
      /* 12 */ coeffs.push_back(mpq_class("1/1100"));
      /* 13 */ coeffs.push_back(mpq_class("1/1188"));
      /* 14 */ coeffs.push_back(mpq_class("1/1276"));
      /* 15 */ coeffs.push_back(mpq_class("1/1364"));
      /* 16 */ coeffs.push_back(mpq_class("1/1452"));
      /* 17 */ coeffs.push_back(mpq_class("1/1540"));
      /* 18 */ coeffs.push_back(mpq_class("1/1628"));
      /* 19 */ coeffs.push_back(mpq_class("1/1716"));
      /* 20 */ coeffs.push_back(mpq_class("1/1804"));
      /* 21 */ coeffs.push_back(mpq_class("1/1892"));
    }
  if (d >= 22)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/46"));
      /* 1  */ coeffs.push_back(mpq_class("1/138"));
      /* 2  */ coeffs.push_back(mpq_class("1/230"));
      /* 3  */ coeffs.push_back(mpq_class("1/322"));
      /* 4  */ coeffs.push_back(mpq_class("1/414"));
      /* 5  */ coeffs.push_back(mpq_class("1/506"));
      /* 6  */ coeffs.push_back(mpq_class("1/598"));
      /* 7  */ coeffs.push_back(mpq_class("1/690"));
      /* 8  */ coeffs.push_back(mpq_class("1/782"));
      /* 9  */ coeffs.push_back(mpq_class("1/874"));
      /* 10 */ coeffs.push_back(mpq_class("1/966"));
      /* 11 */ coeffs.push_back(mpq_class("1/1058"));
      /* 12 */ coeffs.push_back(mpq_class("1/1150"));
      /* 13 */ coeffs.push_back(mpq_class("1/1242"));
      /* 14 */ coeffs.push_back(mpq_class("1/1334"));
      /* 15 */ coeffs.push_back(mpq_class("1/1426"));
      /* 16 */ coeffs.push_back(mpq_class("1/1518"));
      /* 17 */ coeffs.push_back(mpq_class("1/1610"));
      /* 18 */ coeffs.push_back(mpq_class("1/1702"));
      /* 19 */ coeffs.push_back(mpq_class("1/1794"));
      /* 20 */ coeffs.push_back(mpq_class("1/1886"));
      /* 21 */ coeffs.push_back(mpq_class("1/1978"));
      /* 22 */ coeffs.push_back(mpq_class("1/2070"));
    }
  if (d >= 23)
    {
      /* 0  */ coeffs.push_back(mpq_class("1/48"));
      /* 1  */ coeffs.push_back(mpq_class("1/144"));
      /* 2  */ coeffs.push_back(mpq_class("1/240"));
      /* 3  */ coeffs.push_back(mpq_class("1/336"));
      /* 4  */ coeffs.push_back(mpq_class("1/432"));
      /* 5  */ coeffs.push_back(mpq_class("1/528"));
      /* 6  */ coeffs.push_back(mpq_class("1/624"));
      /* 7  */ coeffs.push_back(mpq_class("1/720"));
      /* 8  */ coeffs.push_back(mpq_class("1/816"));
      /* 9  */ coeffs.push_back(mpq_class("1/912"));
      /* 10 */ coeffs.push_back(mpq_class("1/1008"));
      /* 11 */ coeffs.push_back(mpq_class("1/1104"));
      /* 12 */ coeffs.push_back(mpq_class("1/1200"));
      /* 13 */ coeffs.push_back(mpq_class("1/1296"));
      /* 14 */ coeffs.push_back(mpq_class("1/1392"));
      /* 15 */ coeffs.push_back(mpq_class("1/1488"));
      /* 16 */ coeffs.push_back(mpq_class("1/1584"));
      /* 17 */ coeffs.push_back(mpq_class("1/1680"));
      /* 18 */ coeffs.push_back(mpq_class("1/1776"));
      /* 19 */ coeffs.push_back(mpq_class("1/1872"));
      /* 20 */ coeffs.push_back(mpq_class("1/1968"));
      /* 21 */ coeffs.push_back(mpq_class("1/2064"));
      /* 22 */ coeffs.push_back(mpq_class("1/2160"));
      /* 23 */ coeffs.push_back(mpq_class("1/2256"));
    }
  if (d >= 24)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/50"));
      /* 1 */ coeffs.push_back(mpq_class("1/150"));
      /* 2 */ coeffs.push_back(mpq_class("1/250"));
      /* 3 */ coeffs.push_back(mpq_class("1/350"));
      /* 4 */ coeffs.push_back(mpq_class("1/450"));
      /* 5 */ coeffs.push_back(mpq_class("1/550"));
      /* 6 */ coeffs.push_back(mpq_class("1/650"));
      /* 7 */ coeffs.push_back(mpq_class("1/750"));
      /* 8 */ coeffs.push_back(mpq_class("1/850"));
      /* 9 */ coeffs.push_back(mpq_class("1/950"));
      /* 10 */ coeffs.push_back(mpq_class("1/1050"));
      /* 11 */ coeffs.push_back(mpq_class("1/1150"));
      /* 12 */ coeffs.push_back(mpq_class("1/1250"));
      /* 13 */ coeffs.push_back(mpq_class("1/1350"));
      /* 14 */ coeffs.push_back(mpq_class("1/1450"));
      /* 15 */ coeffs.push_back(mpq_class("1/1550"));
      /* 16 */ coeffs.push_back(mpq_class("1/1650"));
      /* 17 */ coeffs.push_back(mpq_class("1/1750"));
      /* 18 */ coeffs.push_back(mpq_class("1/1850"));
      /* 19 */ coeffs.push_back(mpq_class("1/1950"));
      /* 20 */ coeffs.push_back(mpq_class("1/2050"));
      /* 21 */ coeffs.push_back(mpq_class("1/2150"));
      /* 22 */ coeffs.push_back(mpq_class("1/2250"));
      /* 23 */ coeffs.push_back(mpq_class("1/2350"));
      /* 24 */ coeffs.push_back(mpq_class("1/2450"));
    }
  if (d >= 25)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/52"));
      /* 1 */ coeffs.push_back(mpq_class("1/156"));
      /* 2 */ coeffs.push_back(mpq_class("1/260"));
      /* 3 */ coeffs.push_back(mpq_class("1/364"));
      /* 4 */ coeffs.push_back(mpq_class("1/468"));
      /* 5 */ coeffs.push_back(mpq_class("1/572"));
      /* 6 */ coeffs.push_back(mpq_class("1/676"));
      /* 7 */ coeffs.push_back(mpq_class("1/780"));
      /* 8 */ coeffs.push_back(mpq_class("1/884"));
      /* 9 */ coeffs.push_back(mpq_class("1/988"));
      /* 10 */ coeffs.push_back(mpq_class("1/1092"));
      /* 11 */ coeffs.push_back(mpq_class("1/1196"));
      /* 12 */ coeffs.push_back(mpq_class("1/1300"));
      /* 13 */ coeffs.push_back(mpq_class("1/1404"));
      /* 14 */ coeffs.push_back(mpq_class("1/1508"));
      /* 15 */ coeffs.push_back(mpq_class("1/1612"));
      /* 16 */ coeffs.push_back(mpq_class("1/1716"));
      /* 17 */ coeffs.push_back(mpq_class("1/1820"));
      /* 18 */ coeffs.push_back(mpq_class("1/1924"));
      /* 19 */ coeffs.push_back(mpq_class("1/2028"));
      /* 20 */ coeffs.push_back(mpq_class("1/2132"));
      /* 21 */ coeffs.push_back(mpq_class("1/2236"));
      /* 22 */ coeffs.push_back(mpq_class("1/2340"));
      /* 23 */ coeffs.push_back(mpq_class("1/2444"));
      /* 24 */ coeffs.push_back(mpq_class("1/2548"));
      /* 25 */ coeffs.push_back(mpq_class("1/2652"));
    }
  if (d >= 26)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/54"));
      /* 1 */ coeffs.push_back(mpq_class("1/162"));
      /* 2 */ coeffs.push_back(mpq_class("1/270"));
      /* 3 */ coeffs.push_back(mpq_class("1/378"));
      /* 4 */ coeffs.push_back(mpq_class("1/486"));
      /* 5 */ coeffs.push_back(mpq_class("1/594"));
      /* 6 */ coeffs.push_back(mpq_class("1/702"));
      /* 7 */ coeffs.push_back(mpq_class("1/810"));
      /* 8 */ coeffs.push_back(mpq_class("1/918"));
      /* 9 */ coeffs.push_back(mpq_class("1/1026"));
      /* 10 */ coeffs.push_back(mpq_class("1/1134"));
      /* 11 */ coeffs.push_back(mpq_class("1/1242"));
      /* 12 */ coeffs.push_back(mpq_class("1/1350"));
      /* 13 */ coeffs.push_back(mpq_class("1/1458"));
      /* 14 */ coeffs.push_back(mpq_class("1/1566"));
      /* 15 */ coeffs.push_back(mpq_class("1/1674"));
      /* 16 */ coeffs.push_back(mpq_class("1/1782"));
      /* 17 */ coeffs.push_back(mpq_class("1/1890"));
      /* 18 */ coeffs.push_back(mpq_class("1/1998"));
      /* 19 */ coeffs.push_back(mpq_class("1/2106"));
      /* 20 */ coeffs.push_back(mpq_class("1/2214"));
      /* 21 */ coeffs.push_back(mpq_class("1/2322"));
      /* 22 */ coeffs.push_back(mpq_class("1/2430"));
      /* 23 */ coeffs.push_back(mpq_class("1/2538"));
      /* 24 */ coeffs.push_back(mpq_class("1/2646"));
      /* 25 */ coeffs.push_back(mpq_class("1/2754"));
      /* 26 */ coeffs.push_back(mpq_class("1/2862"));
    }
  if (d >= 27)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/56"));
      /* 1 */ coeffs.push_back(mpq_class("1/168"));
      /* 2 */ coeffs.push_back(mpq_class("1/280"));
      /* 3 */ coeffs.push_back(mpq_class("1/392"));
      /* 4 */ coeffs.push_back(mpq_class("1/504"));
      /* 5 */ coeffs.push_back(mpq_class("1/616"));
      /* 6 */ coeffs.push_back(mpq_class("1/728"));
      /* 7 */ coeffs.push_back(mpq_class("1/840"));
      /* 8 */ coeffs.push_back(mpq_class("1/952"));
      /* 9 */ coeffs.push_back(mpq_class("1/1064"));
      /* 10 */ coeffs.push_back(mpq_class("1/1176"));
      /* 11 */ coeffs.push_back(mpq_class("1/1288"));
      /* 12 */ coeffs.push_back(mpq_class("1/1400"));
      /* 13 */ coeffs.push_back(mpq_class("1/1512"));
      /* 14 */ coeffs.push_back(mpq_class("1/1624"));
      /* 15 */ coeffs.push_back(mpq_class("1/1736"));
      /* 16 */ coeffs.push_back(mpq_class("1/1848"));
      /* 17 */ coeffs.push_back(mpq_class("1/1960"));
      /* 18 */ coeffs.push_back(mpq_class("1/2072"));
      /* 19 */ coeffs.push_back(mpq_class("1/2184"));
      /* 20 */ coeffs.push_back(mpq_class("1/2296"));
      /* 21 */ coeffs.push_back(mpq_class("1/2408"));
      /* 22 */ coeffs.push_back(mpq_class("1/2520"));
      /* 23 */ coeffs.push_back(mpq_class("1/2632"));
      /* 24 */ coeffs.push_back(mpq_class("1/2744"));
      /* 25 */ coeffs.push_back(mpq_class("1/2856"));
      /* 26 */ coeffs.push_back(mpq_class("1/2968"));
      /* 27 */ coeffs.push_back(mpq_class("1/3080"));
    }
  if (d >= 28)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/58"));
      /* 1 */ coeffs.push_back(mpq_class("1/174"));
      /* 2 */ coeffs.push_back(mpq_class("1/290"));
      /* 3 */ coeffs.push_back(mpq_class("1/406"));
      /* 4 */ coeffs.push_back(mpq_class("1/522"));
      /* 5 */ coeffs.push_back(mpq_class("1/638"));
      /* 6 */ coeffs.push_back(mpq_class("1/754"));
      /* 7 */ coeffs.push_back(mpq_class("1/870"));
      /* 8 */ coeffs.push_back(mpq_class("1/986"));
      /* 9 */ coeffs.push_back(mpq_class("1/1102"));
      /* 10 */ coeffs.push_back(mpq_class("1/1218"));
      /* 11 */ coeffs.push_back(mpq_class("1/1334"));
      /* 12 */ coeffs.push_back(mpq_class("1/1450"));
      /* 13 */ coeffs.push_back(mpq_class("1/1566"));
      /* 14 */ coeffs.push_back(mpq_class("1/1682"));
      /* 15 */ coeffs.push_back(mpq_class("1/1798"));
      /* 16 */ coeffs.push_back(mpq_class("1/1914"));
      /* 17 */ coeffs.push_back(mpq_class("1/2030"));
      /* 18 */ coeffs.push_back(mpq_class("1/2146"));
      /* 19 */ coeffs.push_back(mpq_class("1/2262"));
      /* 20 */ coeffs.push_back(mpq_class("1/2378"));
      /* 21 */ coeffs.push_back(mpq_class("1/2494"));
      /* 22 */ coeffs.push_back(mpq_class("1/2610"));
      /* 23 */ coeffs.push_back(mpq_class("1/2726"));
      /* 24 */ coeffs.push_back(mpq_class("1/2842"));
      /* 25 */ coeffs.push_back(mpq_class("1/2958"));
      /* 26 */ coeffs.push_back(mpq_class("1/3074"));
      /* 27 */ coeffs.push_back(mpq_class("1/3190"));
      /* 28 */ coeffs.push_back(mpq_class("1/3306"));
    }
  if (d >= 29)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/60"));
      /* 1 */ coeffs.push_back(mpq_class("1/180"));
      /* 2 */ coeffs.push_back(mpq_class("1/300"));
      /* 3 */ coeffs.push_back(mpq_class("1/420"));
      /* 4 */ coeffs.push_back(mpq_class("1/540"));
      /* 5 */ coeffs.push_back(mpq_class("1/660"));
      /* 6 */ coeffs.push_back(mpq_class("1/780"));
      /* 7 */ coeffs.push_back(mpq_class("1/900"));
      /* 8 */ coeffs.push_back(mpq_class("1/1020"));
      /* 9 */ coeffs.push_back(mpq_class("1/1140"));
      /* 10 */ coeffs.push_back(mpq_class("1/1260"));
      /* 11 */ coeffs.push_back(mpq_class("1/1380"));
      /* 12 */ coeffs.push_back(mpq_class("1/1500"));
      /* 13 */ coeffs.push_back(mpq_class("1/1620"));
      /* 14 */ coeffs.push_back(mpq_class("1/1740"));
      /* 15 */ coeffs.push_back(mpq_class("1/1860"));
      /* 16 */ coeffs.push_back(mpq_class("1/1980"));
      /* 17 */ coeffs.push_back(mpq_class("1/2100"));
      /* 18 */ coeffs.push_back(mpq_class("1/2220"));
      /* 19 */ coeffs.push_back(mpq_class("1/2340"));
      /* 20 */ coeffs.push_back(mpq_class("1/2460"));
      /* 21 */ coeffs.push_back(mpq_class("1/2580"));
      /* 22 */ coeffs.push_back(mpq_class("1/2700"));
      /* 23 */ coeffs.push_back(mpq_class("1/2820"));
      /* 24 */ coeffs.push_back(mpq_class("1/2940"));
      /* 25 */ coeffs.push_back(mpq_class("1/3060"));
      /* 26 */ coeffs.push_back(mpq_class("1/3180"));
      /* 27 */ coeffs.push_back(mpq_class("1/3300"));
      /* 28 */ coeffs.push_back(mpq_class("1/3420"));
      /* 29 */ coeffs.push_back(mpq_class("1/3540"));
    }
  if (d >= 30)
    {
      /* 0 */ coeffs.push_back(mpq_class("1/62"));
      /* 1 */ coeffs.push_back(mpq_class("1/186"));
      /* 2 */ coeffs.push_back(mpq_class("1/310"));
      /* 3 */ coeffs.push_back(mpq_class("1/434"));
      /* 4 */ coeffs.push_back(mpq_class("1/558"));
      /* 5 */ coeffs.push_back(mpq_class("1/682"));
      /* 6 */ coeffs.push_back(mpq_class("1/806"));
      /* 7 */ coeffs.push_back(mpq_class("1/930"));
      /* 8 */ coeffs.push_back(mpq_class("1/1054"));
      /* 9 */ coeffs.push_back(mpq_class("1/1178"));
      /* 10 */ coeffs.push_back(mpq_class("1/1302"));
      /* 11 */ coeffs.push_back(mpq_class("1/1426"));
      /* 12 */ coeffs.push_back(mpq_class("1/1550"));
      /* 13 */ coeffs.push_back(mpq_class("1/1674"));
      /* 14 */ coeffs.push_back(mpq_class("1/1798"));
      /* 15 */ coeffs.push_back(mpq_class("1/1922"));
      /* 16 */ coeffs.push_back(mpq_class("1/2046"));
      /* 17 */ coeffs.push_back(mpq_class("1/2170"));
      /* 18 */ coeffs.push_back(mpq_class("1/2294"));
      /* 19 */ coeffs.push_back(mpq_class("1/2418"));
      /* 20 */ coeffs.push_back(mpq_class("1/2542"));
      /* 21 */ coeffs.push_back(mpq_class("1/2666"));
      /* 22 */ coeffs.push_back(mpq_class("1/2790"));
      /* 23 */ coeffs.push_back(mpq_class("1/2914"));
      /* 24 */ coeffs.push_back(mpq_class("1/3038"));
      /* 25 */ coeffs.push_back(mpq_class("1/3162"));
      /* 26 */ coeffs.push_back(mpq_class("1/3286"));
      /* 27 */ coeffs.push_back(mpq_class("1/3410"));
      /* 28 */ coeffs.push_back(mpq_class("1/3534"));
      /* 29 */ coeffs.push_back(mpq_class("1/3658"));
      /* 30 */ coeffs.push_back(mpq_class("1/3782"));
    }
}
