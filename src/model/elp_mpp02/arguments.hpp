// Parameter adjustments and fundamental arguments adapted from Yuk Tung Liu's
// ELP/MPP02 implementation, revision e25e08f6082a3f3211b3d14e062d935558ad3ea6.
// https://github.com/ytliu0/ElpMpp02 (GPL-3.0). See data/elp_mpp02/README.md.
// Changes: internal names/namespace, standard math names and immutable callers.
#pragma once

#include <cmath>
#include "spp/core/constants.hpp"

namespace spp::model::elp_mpp02::detail {

inline double mod2pi(double x) {
    return x - two_pi * std::floor((x + pi) / two_pi);
}

struct Arguments {
  double W1, D, F, L, Lp, zeta, Me, Ve, EM, Ma, Ju, Sa, Ur, Ne;
};

// Factors multiplied by B1-B5 for longitude and latitude
struct Factors {
  double fA, fB1, fB2, fB3, fB4, fB5;
};

// parameters for adjusting the lunar and planetary arguments
struct Parameters {
  // parameters adjusted to fit data
  double Dw1_0, Dw2_0, Dw3_0, Deart_0, Dperi, Dw1_1, Dgam, De, Deart_1, Dep,
         Dw2_1, Dw3_1, Dw1_2, Dw1_3, Dw1_4, Dw2_2, Dw2_3, Dw3_2, Dw3_3;
  // parameters derived from the previous parameters
  double Cw2_1, Cw3_1;
};

void setup_parameters(int corr, Parameters &paras, Factors &facs) {
   // PARAMETERS adjusted to fit data
   switch (corr) {
       case 0:
         paras.Dw1_0   = -0.10525;
         paras.Dw2_0   =  0.16826;
         paras.Dw3_0   = -0.10760;
         paras.Deart_0 = -0.04012;
         paras.Dperi   = -0.04854;
         paras.Dw1_1   = -0.32311;
         paras.Dgam    =  0.00069;
         paras.De      =  0.00005;
         paras.Deart_1 =  0.01442;
         paras.Dep     =  0.00226;
         paras.Dw2_1   =  0.08017;
         paras.Dw3_1   = -0.04317;
         paras.Dw1_2   = -0.03794;
         paras.Dw1_3   =  0.0;
         paras.Dw1_4   =  0.0;
         paras.Dw2_2   =  0.0;
         paras.Dw2_3   =  0.0;
         paras.Dw3_2   =  0.0;
         paras.Dw3_3   =  0.0;
         break;
      case 1:
         paras.Dw1_0   = -0.07008;
         paras.Dw2_0   =  0.20794;
         paras.Dw3_0   = -0.07215;
         paras.Deart_0 = -0.00033;
         paras.Dperi   = -0.00749;
         paras.Dw1_1   = -0.35106;
         paras.Dgam    =  0.00085;
         paras.De      = -0.00006;
         paras.Deart_1 =  0.00732;
         paras.Dep     =  0.00224;
         paras.Dw2_1   =  0.08017;
         paras.Dw3_1   = -0.04317;
         paras.Dw1_2   = -0.03743;
         paras.Dw1_3   = -0.00018865;
         paras.Dw1_4   = -0.00001024;
         paras.Dw2_2   =  0.00470602;
         paras.Dw2_3   = -0.00025213;
         paras.Dw3_2   = -0.00261070;
         paras.Dw3_3   = -0.00010712;
         break;
   }
   // derived parameters
   const double am = 0.074801329;
   const double alpha = 0.002571881;
   const double dtsm = 2.0*alpha/(3.0*am);
   const double xa = 2.0*alpha/3.0;
   const double sec = pi/648000.0; // arcsecs -> radians
   double bp[5][2] = {{0.311079095, -0.103837907},
                      {-0.004482398, 0.000668287},
                      {-0.001102485, -0.001298072},
                      {0.001056062, -0.000178028},
                      {0.000050928, -0.000037342}};
   double w11 = (1732559343.73604 + paras.Dw1_1)*sec;
   double w21 = (14643420.3171 + paras.Dw2_1)*sec;
   double w31 = (-6967919.5383 + paras.Dw3_1)*sec;
   double x2 = w21/w11;
   double x3 = w31/w11;
   double y2 = am*bp[0][0] + xa*bp[4][0];
   double y3 = am*bp[0][1] + xa*bp[4][1];
   double d21 = x2-y2;
   double d22 = w11*bp[1][0];
   double d23 = w11*bp[2][0];
   double d24 = w11*bp[3][0];
   double d25 = y2/am;
   double d31 = x3-y3;
   double d32 = w11*bp[1][1];
   double d33 = w11*bp[2][1];
   double d34 = w11*bp[3][1];
   double d35 = y3/am;
   paras.Cw2_1 = d21*paras.Dw1_1 + d25*paras.Deart_1 + d22*paras.Dgam +
                  d23*paras.De + d24*paras.Dep;
   paras.Cw3_1 = d31*paras.Dw1_1 + d35*paras.Deart_1 + d32*paras.Dgam +
                  d33*paras.De + d34*paras.Dep;

   // factors multipled by B1-B5 for longitude and latitude
   double delnu_nu = (0.55604 + paras.Dw1_1)*sec/w11;
   double dele = (0.01789 + paras.De)*sec;
   double delg = (-0.08066 + paras.Dgam)*sec;
   double delnp_nu = (-0.06424 + paras.Deart_1)*sec/w11;
   double delep = (-0.12879 + paras.Dep)*sec;
   // factors multipled by B1-B5 for longitude and latitude
   facs.fB1 = -am*delnu_nu + delnp_nu;
   facs.fB2 = delg;
   facs.fB3 = dele;
   facs.fB4 = delep;
   facs.fB5 = -xa*delnu_nu + dtsm*delnp_nu;
   // factor multiplie A_i for distance
   facs.fA = 1.0 - 2.0/3.0*delnu_nu;
}

void compute_Elp_arguments(double T, Parameters paras, Arguments &args) {
     const double deg = pi/180.0; // degrees -> radians
     const double sec = pi/648000.0; // arcsecs -> radians
     double T2 = T*T;
     double T3 = T*T2;
     double T4 = T2*T2;
     double w10 = (-142.0 + 18.0/60.0 +(59.95571 + paras.Dw1_0)/3600.0)*deg;
     double w11 = mod2pi((1732559343.73604 + paras.Dw1_1)*T*sec);
     double w12 = mod2pi((-6.8084 + paras.Dw1_2)*T2*sec);
     double w13 = mod2pi((0.006604 + paras.Dw1_3)*T3*sec);
     double w14 = mod2pi((-3.169e-5 + paras.Dw1_4)*T4*sec);
     double w20 = (83.0 + 21.0/60.0 + (11.67475 + paras.Dw2_0)/3600.0)*deg;
     double w21 = mod2pi((14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*T*sec);
     double w22 = mod2pi((-38.2631 + paras.Dw2_2)*T2*sec);
     double w23 = mod2pi((-0.045047+ paras.Dw2_3)*T3*sec);
     double w24 = mod2pi(0.00021301*T4*sec);
     double w30 = (125.0 + 2.0/60.0 + (40.39816 + paras.Dw3_0)/3600.0)*deg;
     double w31 = mod2pi((-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*T*sec);
     double w32 = mod2pi((6.359 + paras.Dw3_2)*T2*sec);
     double w33 = mod2pi((0.007625 + paras.Dw3_3)*T3*sec);
     double w34 = mod2pi(-3.586e-5*T4*sec);
     double Ea0 = (100.0 + 27.0/60.0 + (59.13885 + paras.Deart_0)/3600.0)*deg;
     double Ea1 = mod2pi((129597742.293 + paras.Deart_1)*T*sec);
     double Ea2 = mod2pi(-0.0202*T2*sec);
     double Ea3 = mod2pi(9e-6*T3*sec);
     double Ea4 = mod2pi(1.5e-7*T4*sec);
     double p0 = (102.0 + 56.0/60.0 + (14.45766 + paras.Dperi)/3600.0)*deg;
     double p1 = mod2pi(1161.24342*T*sec);
     double p2 = mod2pi(0.529265*T2*sec);
     double p3 = mod2pi(-1.1814e-4*T3*sec);
     double p4 = mod2pi(1.1379e-5*T4*sec);

     double Me = (-108.0 + 15.0/60.0 + 3.216919/3600.0)*deg;
     Me += mod2pi(538101628.66888*T*sec);
     double Ve = (-179.0 + 58.0/60.0 + 44.758419/3600.0)*deg;
     Ve += mod2pi(210664136.45777*T*sec);
     double EM = (100.0 + 27.0/60.0 + 59.13885/3600.0)*deg;
     EM += mod2pi(129597742.293*T*sec);
     double Ma = (-5.0 + 26.0/60.0 + 3.642778/3600.0)*deg;
     Ma += mod2pi(68905077.65936*T*sec);
     double Ju = (34.0 + 21.0/60.0 + 5.379392/3600.0)*deg;
     Ju += mod2pi(10925660.57335*T*sec);
     double Sa = (50.0 + 4.0/60.0 + 38.902495/3600.0)*deg;
     Sa += mod2pi(4399609.33632*T*sec);
     double Ur = (-46.0 + 3.0/60.0 + 4.354234/3600.0)*deg;
     Ur += mod2pi(1542482.57845*T*sec);
     double Ne = (-56.0 + 20.0/60.0 + 56.808371/3600.0)*deg;
     Ne += mod2pi(786547.897*T*sec);

     double W1 = w10+w11+w12+w13+w14;
     double W2 = w20+w21+w22+w23+w24;
     double W3 = w30+w31+w32+w33+w34;
     double Ea = Ea0+Ea1+Ea2+Ea3+Ea4;
     double pomp = p0+p1+p2+p3+p4;

     // Mean longitude of the Moon
     args.W1 = mod2pi(W1);
     // Arguments of Delaunay
     args.D = mod2pi(W1-Ea + pi);
     args.F = mod2pi(W1-W3);
     args.L = mod2pi(W1-W2);
     args.Lp = mod2pi(Ea-pomp);

     // zeta
     args.zeta = mod2pi(W1 + 0.02438029560881907*T);

     // Planetary arguments (mean longitudes and mean motions)
     args.Me = mod2pi(Me);
     args.Ve = mod2pi(Ve);
     args.EM = mod2pi(EM);
     args.Ma = mod2pi(Ma);
     args.Ju = mod2pi(Ju);
     args.Sa = mod2pi(Sa);
     args.Ur = mod2pi(Ur);
     args.Ne = mod2pi(Ne);
}

// Compute the time derivatives of the lunar and planetary arguments used in the ELP/MPP02 series
// return the time derivatives in units of rad/day
void compute_Elp_arguments_dot(double T, Parameters paras, Arguments &args_dot) {
     const double fac = pi/648000.0/36525; // arcsecs -> radians/cy
     double T2 = T*T;
     double T3 = T*T2;
     double w11 = (1732559343.73604 + paras.Dw1_1)*fac;
     double w12 = (-6.8084 + paras.Dw1_2)*2*T*fac;
     double w13 = (0.006604 + paras.Dw1_3)*3*T2*fac;
     double w14 = (-3.169e-5 + paras.Dw1_4)*4*T3*fac;
     double w21 = (14643420.3171 + paras.Dw2_1 + paras.Cw2_1)*fac;
     double w22 = (-38.2631 + paras.Dw2_2)*2*T*fac;
     double w23 = (-0.045047+ paras.Dw2_3)*3*T2*fac;
     double w24 = 0.00021301*4*T3*fac;
     double w31 = (-6967919.5383 + paras.Dw3_1 + paras.Cw3_1)*fac;
     double w32 = (6.359 + paras.Dw3_2)*2*T*fac;
     double w33 = (0.007625 + paras.Dw3_3)*3*T2*fac;
     double w34 = -3.586e-5*4*T3*fac;
     double Ea1 = (129597742.293 + paras.Deart_1)*fac;
     double Ea2 = -0.0202*2*T*fac;
     double Ea3 = 9e-6*3*T2*fac;
     double Ea4 = 1.5e-7*4*T3*fac;
     double p1 = 1161.24342*fac;
     double p2 = 0.529265*2*T*fac;
     double p3 = -1.1814e-4*3*T2*fac;
     double p4 = 1.1379e-5*4*T3*fac;

     double W1 = w11+w12+w13+w14;
     double W2 = w21+w22+w23+w24;
     double W3 = w31+w32+w33+w34;
     double Ea = Ea1+Ea2+Ea3+Ea4;
     double pomp = p1+p2+p3+p4;

     // Mean longitude of the Moon
     args_dot.W1 = W1;
     // Arguments of Delaunay
     args_dot.D = W1-Ea;
     args_dot.F = W1-W3;
     args_dot.L = W1-W2;
     args_dot.Lp = Ea-pomp;

     // zeta
     args_dot.zeta = W1 + 5028.79695*fac;

     // Planetary arguments (mean longitudes and mean motions)
     args_dot.Me = 538101628.66888*fac;
     args_dot.Ve = 210664136.45777*fac;
     args_dot.EM = 129597742.293*fac;
     args_dot.Ma = 68905077.65936*fac;
     args_dot.Ju = 10925660.57335*fac;
     args_dot.Sa = 4399609.33632*fac;
     args_dot.Ur = 1542482.57845*fac;
     args_dot.Ne = 786547.897*fac;
}


}  // namespace spp::model::elp_mpp02::detail
