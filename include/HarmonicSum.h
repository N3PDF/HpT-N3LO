/*
 * =====================================================================================
 *
 *       Filename:
 *
 *    Description: 
 *
 *        Version:  1.0
 *        Created:  
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */


#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <sstream>
#include <string>
#include <functional>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_zeta.h>

#include "ComplexDefs.h"

using namespace std;

class HSum
{
    public:
        HSum(bool verbose=false,bool testinterfun=false, bool testharmsums=false);
        virtual ~HSum();

        std::complex<double> HS(int i, std::complex<double> N);
        std::complex<double> HS(int i, int j, std::complex<double> N);
        std::complex<double> HS(int i, int j, int k, std::complex<double> N);
        std::complex<double> HS(int i, int j, int k, int m, std::complex<double> N);

    private:
        bool _testinterpolatedfunction;
        bool _testharmonicsums;
        bool _verbose;

        void InizializeConst();

        // crosscheck function
        void TestInterpolatedFunction();
        void TestHarmonicSums(int n);
        std::complex<double> HS_int(int i, int N);
        std::complex<double> HS_int(int i, int j, int N);
        std::complex<double> HS_int(int i, int j, int k, int N);
        std::complex<double> HS_int(int i, int j, int k, int m, int N);

        // Semi-analitics Mellin transform
        std::complex <double> g1(std::complex<double> N);
        std::complex <double> g2(std::complex<double> N);
        std::complex <double> g3(std::complex<double> N);
        std::complex <double> g4(std::complex<double> N);
        std::complex <double> g5(std::complex<double> N);
        std::complex <double> g6(std::complex<double> N);
        std::complex <double> g7(std::complex<double> N);
        std::complex <double> g8(std::complex<double> N);
        std::complex <double> g9(std::complex<double> N);
        std::complex <double> g10(std::complex<double> N);
        std::complex <double> g11(std::complex<double> N);
        std::complex <double> g12(std::complex<double> N);
        std::complex <double> g13(std::complex<double> N);
        std::complex <double> g14(std::complex<double> N);
        std::complex <double> g15(std::complex<double> N);
        std::complex <double> g16(std::complex<double> N);
        std::complex <double> g17(std::complex<double> N);
        std::complex <double> g18(std::complex<double> N);
        std::complex <double> g19(std::complex<double> N);
        std::complex <double> g20(std::complex<double> N);
        std::complex <double> g21(std::complex<double> N);
        std::complex <double> g22(std::complex<double> N);
        std::complex <double> g23(std::complex<double> N);
        std::complex <double> g24(std::complex<double> N);
        std::complex <double> g25(std::complex<double> N);
        std::complex <double> g26(std::complex<double> N);
        std::complex <double> g27(std::complex<double> N);
        std::complex <double> g28(std::complex<double> N);
        std::complex <double> g29(std::complex<double> N);
        std::complex <double> g30(std::complex<double> N);
        std::complex <double> g31(std::complex<double> N);
        std::complex <double> g32(std::complex<double> N);
        std::complex <double> g33(std::complex<double> N);
        std::complex <double> g34(std::complex<double> N);
        std::complex <double> g35(std::complex<double> N);
        std::complex <double> g36(std::complex<double> N);
        std::complex <double> g37(std::complex<double> N);
        std::complex <double> g38(std::complex<double> N);
        std::complex <double> g39(std::complex<double> N);

        //Math important Function
        double Zeta(int i)
        {
          return gsl_sf_zeta_int(i);
        }

        std::complex<double> B(int i,std::complex<double> z)
        {
          return (1./std::pow(2.,(double)i+1.)*(PolyGamma(i,(1.+z)/2.)
                      -PolyGamma(i,z/2.)));
        }
        //Important Costants
        double zeta2;
        double zeta3;
        double Li4;
        double log2;
        double log2q;
        double log2c;
        double zeta2q;
        double EulerGamma;

        //Vectors of coefficients
        std::vector<double> a1;
        std::vector<double> a2;
        std::vector<double> a3;
        std::vector<double> b1;
        std::vector<double> b2;
        std::vector<double> b3;
        std::vector<double> c1;
        std::vector<double> P21;
        std::vector<double> c2;
        std::vector<double> P22;
        std::vector<double> c3;
        std::vector<double> P23;
        std::vector<double> P33;
        std::vector<double> c4;
        std::vector<double> P24;
        std::vector<double> P34;
        std::vector<double> c5;
        std::vector<double> d5;
        std::vector<double> q1;
        std::vector<double> q2;
        std::vector<double> q3;
        std::vector<double> q4;
        std::vector<double> q5;
        std::vector<double> q6;
        std::vector<double> q7;

        // Table of Harmonic Sums up to weight 4
        // Weight 1
        std::complex<double> H_1(std::complex<double> N);
        std::complex<double> H_m1(std::complex<double> N);
        // Weight 2
        std::complex<double> H_2(std::complex<double> N);
        std::complex<double> H_m2(std::complex<double> N);
        std::complex<double> H_1_1(std::complex<double> N);
        std::complex<double> H_1_m1(std::complex<double> N);
        std::complex<double> H_m1_1(std::complex<double> N);
        std::complex<double> H_m1_m1(std::complex<double> N);
        // weight 3
        std::complex<double> H_3(std::complex<double> N);
        std::complex<double> H_m3(std::complex<double> N);
        std::complex<double> H_1_2(std::complex<double> N);
        std::complex<double> H_2_1(std::complex<double> N);
        std::complex<double> H_m2_1(std::complex<double> N);
        std::complex<double> H_m1_2(std::complex<double> N);
        std::complex<double> H_1_m2(std::complex<double> N);
        std::complex<double> H_2_m1(std::complex<double> N);
        std::complex<double> H_m2_m1(std::complex<double> N);
        std::complex<double> H_m1_m2(std::complex<double> N);
        std::complex<double> H_1_1_1(std::complex<double> N);
        std::complex<double> H_m1_1_1(std::complex<double> N);
        std::complex<double> H_1_1_m1(std::complex<double> N);
        std::complex<double> H_1_m1_1(std::complex<double> N);
        std::complex<double> H_1_m1_m1(std::complex<double> N);
        std::complex<double> H_m1_1_m1(std::complex<double> N);
        std::complex<double> H_m1_m1_1(std::complex<double> N);
        std::complex<double> H_m1_m1_m1(std::complex<double> N);
        // Weight 4
        std::complex<double> H_4(std::complex<double> N);
        std::complex<double> H_m4(std::complex<double> N);
        std::complex<double> H_m3_m1(std::complex<double> N);
        std::complex<double> H_m3_1(std::complex<double> N);
        std::complex<double> H_m1_m3(std::complex<double> N);
        std::complex<double> H_m1_3(std::complex<double> N);
        std::complex<double> H_1_m3(std::complex<double> N);
        std::complex<double> H_1_3(std::complex<double> N);
        std::complex<double> H_3_1(std::complex<double> N);
        std::complex<double> H_3_m1(std::complex<double> N);
        std::complex<double> H_m2_m2(std::complex<double> N);
        std::complex<double> H_m2_2(std::complex<double> N);
        std::complex<double> H_2_m2(std::complex<double> N);
        std::complex<double> H_2_2(std::complex<double> N);
        std::complex<double> H_m2_1_1(std::complex<double> N);
        std::complex<double> H_m2_1_m1(std::complex<double> N);
        std::complex<double> H_m2_m1_1(std::complex<double> N);
        std::complex<double> H_m2_m1_m1(std::complex<double> N);
        std::complex<double> H_m1_1_2(std::complex<double> N);
        std::complex<double> H_m1_1_m2(std::complex<double> N);
        std::complex<double> H_m1_m1_2(std::complex<double> N);
        std::complex<double> H_m1_m1_m2(std::complex<double> N);
        std::complex<double> H_m1_2_1(std::complex<double> N);
        std::complex<double> H_m1_2_m1(std::complex<double> N);
        std::complex<double> H_m1_m2_1(std::complex<double> N);
        std::complex<double> H_m1_m2_m1(std::complex<double> N);
        std::complex<double> H_1_1_2(std::complex<double> N);
        std::complex<double> H_1_1_m2(std::complex<double> N);
        std::complex<double> H_1_m1_2(std::complex<double> N);
        std::complex<double> H_1_m1_m2(std::complex<double> N);
        std::complex<double> H_1_2_1(std::complex<double> N);
        std::complex<double> H_1_2_m1(std::complex<double> N);
        std::complex<double> H_1_m2_1(std::complex<double> N);
        std::complex<double> H_1_m2_m1(std::complex<double> N);
        std::complex<double> H_2_1_1(std::complex<double> N);
        std::complex<double> H_2_1_m1(std::complex<double> N);
        std::complex<double> H_2_m1_1(std::complex<double> N);
        std::complex<double> H_2_m1_m1(std::complex<double> N);
        std::complex<double> H_1_1_1_1(std::complex<double> N);
        std::complex<double> H_1_1_1_m1(std::complex<double> N);
        std::complex<double> H_1_1_m1_1(std::complex<double> N);
        std::complex<double> H_1_1_m1_m1(std::complex<double> N);
        std::complex<double> H_1_m1_1_1(std::complex<double> N);
        std::complex<double> H_1_m1_1_m1(std::complex<double> N);
        std::complex<double> H_1_m1_m1_1(std::complex<double> N);
        std::complex<double> H_1_m1_m1_m1(std::complex<double> N);
        std::complex<double> H_m1_1_1_1(std::complex<double> N);
        std::complex<double> H_m1_1_1_m1(std::complex<double> N);
        std::complex<double> H_m1_1_m1_1(std::complex<double> N);
        std::complex<double> H_m1_1_m1_m1(std::complex<double> N);
        std::complex<double> H_m1_m1_1_1(std::complex<double> N);
        std::complex<double> H_m1_m1_1_m1(std::complex<double> N);
        std::complex<double> H_m1_m1_m1_1(std::complex<double> N);
        std::complex<double> H_m1_m1_m1_m1(std::complex<double> N);

        std::complex<double> Log(std::complex<double> N)
        {
           return(std::log(N));
        }

        //Complex Euler Gamma
        std::complex<double> LogGamma(std::complex<double> z)
        {
          int g=7;
          double p[9];
          p[0]=0.99999999999980993, p[1]=676.5203681218851;
          p[2]=-1259.1392167224028, p[3]=771.32342877765313;
          p[4]=-176.61502916214059, p[5]=12.507343278686905;
          p[6]=-0.1385710952657201, p[7]=9.9843695780195716e-6;
          p[8]=1.5056327351493116e-7;

          std::complex<double>sum,ris;
          z  -= 1;
          sum = p[0];
          for(int i=1;i<(g+2);i++)
            sum += p[i]/(z+i);
          ris=0.5*gsl_sf_log(2*M_PI)+(z+0.5)*log(z+g+0.5)-(z+g+0.5)+log(sum);
          return ris;
        }

        std::complex<double> CGamma(std::complex<double> z)
        {
          return( std::exp(LogGamma(z)));
        }

        //Gamma Derivatives of order i complex
        std::complex<double> PolyGamma(int i, std::complex<double> z)
        {
          if (i == 0)
          {
            std::complex<double> SUB = 0. ;
            std::complex<double> ZZ = z;
            if(std::abs(std::imag(ZZ))<10.)
            { // if too close to the real axis...
            label1:
              if(std::real(ZZ)<10.)
              { // ...use recurrence relation to push real(z) large enough
                SUB = SUB - 1./ ZZ;
                ZZ = ZZ + 1.;
                goto label1;
              }
            }
            std::complex<double> RZ = 1./ ZZ;
            std::complex<double> DZ = RZ * RZ;
            // SUB + asympt expansion (Abramowitz, Stengun, 6.3.18)
            return (SUB+std::log(ZZ)-0.5*RZ-DZ/5040.*( 420.+DZ*(-42.+DZ*(20.-21.*DZ))));
          }
          else
          {
             int K1, K2;
             std::complex<double> SUB = 0. , SUBM;
             std::complex<double> ZZ = z;
             if(std::abs(std::imag(ZZ))<10.)
             { // if too close to the real axis...
               label2:
               SUBM = -1./ZZ;
               for(K1=1; K1<=i; K1++)
               {
                 SUBM = - SUBM * K1 / ZZ;
               }
              if(std::real(ZZ)<10.)
              { // ...use recurrence relation to push real(z) large enough
                SUB = SUB + SUBM;
                ZZ = ZZ + 1.;
                goto label2;
              }
            }

            // Expansion coefficients for the first derivative
            double A1 =  1.;
            double A2 =  1./2.;
            double A3 =  1./6.;
            double A4 = -1./30.;
            double A5 =  1./42.;
            double A6 = -1./30.;
            double A7 =  5./66.;

            // Expansion coefficients for the higher derivatives
            if(i>1)
            {
              for(K2=2;K2<=i;K2++)
              {
                A1 = A1 * (1.*K2-1.);
                A2 = A2 *  1.*K2;
                A3 = A3 * (1.*K2+1.);
                A4 = A4 * (1.*K2+3.);
                A5 = A5 * (1.*K2+5.);
                A6 = A6 * (1.*K2+7.);
                A7 = A7 * (1.*K2+9.);
              }
            }
            std::complex<double> RZ = 1./ ZZ;
            std::complex<double> DZ = RZ * RZ;

            // SUB + asympt expansion (Abramowitz, Stengun, 6.4.11)
            return (SUB+std::pow(-1.,i+1)*std::pow(RZ,i)*
                    (A1+RZ*(A2+RZ*(A3+DZ*(A4+DZ*(A5+DZ*(A6+A7*DZ)))))));
          }
        }
};
