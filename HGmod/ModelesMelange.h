// Copyright (c) 2016 Bernard Giroux. All rights reserved.
/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef __MODELESMELANGE_H__
#define __MODELESMELANGE_H__

/*
 * References
 *
 *
 * @inproceedings{druchinin00,
 *   author = {S. V. Druchinin},
 *   address = {Gold Coast, Australia},
 *   bootitle = {Proceedings of the 8$^{th}$ International Conference on Ground-Penetrating Radar},
 *   editor = {D. Noon and G. Stickley and D. Longstaff},
 *   title = {Models for calculation of dielectric constant of moist sandy-clayey soils in wavelengths from centimeters to tens of meters},
 *   year = {2000}
 *}
 *
 * @article{sihvola00,
 *   author = {Ari Sihvola},
 *   journal = ssta,
 *   number = {4},
 *   pages = {393--415},
 *   title = {Mixing Rules with Complex Dielectric Coefficients},
 *   volume = {1},
 *   year = {2000}
 * }
 *
 * @book{hasted73,
 *   author = {J. B. Hasted},
 *   address = {London},
 *   publisher = {Chapman and Hall},
 *   title = {Aqueous Dielectrics},
 *   year = {1973}
 * }
 *
 * @INCOLLECTION{olhoeft81,
 *     author = {Gary R. Olhoeft},
 *     title = {Electrical properties of rocks},
 *     booktitle = {Physical properties of Rocks and Minerals},
 *     publisher = {McGraw-Hill},
 *     year = {1981},
 *     editor = {Y. S. Touloukian and Y. S. Judd and R. F. Roy},
 *     pages = {257--330},
 *     address = {New York},
 * }
 *
 *
 * @Article{chang11,
 *   Title   = {A Parallel Derivation to the Maxwell-Garnett Formula for the Magnetic Permeability of Mixed Materials},
 *   Author  = {Hsien-Ming Chang and Chungpin Liao},
 *   Journal = {World Journal of Condensed Matter Physics},
 *   Year    = {2011},
 *   Number  = {2},
 *   Pages   = {55--58},
 *   Volume  = {1},
 *   Doi     = {10.4236/wjcmp.2011.12009}
 * }
 *
 */

#include <cmath>
#include <complex>
#include <valarray>
#include <vector>
#include <iostream>
#include <algorithm>

#include "constantes.h"

#include "Powell.h"

namespace ModelesMelange
{
    template<typename T>
    T sign(const T x) { return x>=0. ? 1. : -1.; }
    
    // eau
    //   const double kappa_0_w = 80.1;
    //   const double kappa_i_w = 4.23;
    //   const double tau_w = 9.3e-12;
    //   const double f_w = 1./tau_w;
    //   const double q_w = 0.987;
    
    // permittivite statique de l'eau en fct de la temperature
    // eq 2.8 de Hasted (1973)
    template<typename T>
    T MalmbergMaryott(const T t) { // t -> temperature en C
        return 87.74 - 0.40008*t + 9.398e-4*t*t - 1.41e-6*t*t*t;
    }
    
    // tableau 2.2 de Hasted (1973)
    template<typename T>
    T relax_w_epsilon_i(const T t) { // t -> temperature en C
        T d[] = {4.46, 4.10, 4.23, 4.20, 4.16, 4.13, 4.21, 4.49};
        if ( t<=0. )
            return d[0];
        else if ( t<10. )
            return d[0] + 0.1*(d[1]-d[0])*t;
        else if ( t<20. )
            return d[1] + 0.1*(d[2]-d[1])*(t-10.);
        else if ( t<30. )
            return d[2] + 0.1*(d[3]-d[2])*(t-20.);
        else if ( t<40. )
            return d[3] + 0.1*(d[4]-d[3])*(t-30.);
        else if ( t<50. )
            return d[4] + 0.1*(d[5]-d[4])*(t-40.);
        else if ( t<60. )
            return d[5] + 0.1*(d[6]-d[5])*(t-50.);
        else if ( t<75. )
            return d[6] + 0.06666666*(d[7]-d[6])*(t-60.);
        else
            return d[7];
    }
    
    // tableau 2.2 de Hasted (1973)
    template<typename T>
    T relax_w_tau(const T t) { // t -> temperature en C
        T d[] = {1.79e-11, 1.26e-11, 0.93e-11, 0.72e-11, 0.58e-11, 0.48e-11,
            0.39e-11, 0.32e-11};
        if ( t<=0. )
            return d[0];
        else if ( t<10. )
            return d[0] + 0.1*(d[1]-d[0])*t;
        else if ( t<20. )
            return d[1] + 0.1*(d[2]-d[1])*(t-10.);
        else if ( t<30. )
            return d[2] + 0.1*(d[3]-d[2])*(t-20.);
        else if ( t<40. )
            return d[3] + 0.1*(d[4]-d[3])*(t-30.);
        else if ( t<50. )
            return d[4] + 0.1*(d[5]-d[4])*(t-40.);
        else if ( t<60. )
            return d[5] + 0.1*(d[6]-d[5])*(t-50.);
        else if ( t<75. )
            return d[6] + 0.06666666*(d[7]-d[6])*(t-60.);
        else
            return d[7];
    }
    
    // tableau 2.2 de Hasted (1973)
    template<typename T>
    T relax_w_alpha(const T t) { // t -> temperature en C
        T d[] = {0.014, 0.014, 0.013, 0.012, 0.009, 0.013, 0.011};
        if ( t<=0. )
            return d[0];
        else if ( t<10. )
            return d[0] + 0.1*(d[1]-d[0])*t;
        else if ( t<20. )
            return d[1] + 0.1*(d[2]-d[1])*(t-10.);
        else if ( t<30. )
            return d[2] + 0.1*(d[3]-d[2])*(t-20.);
        else if ( t<40. )
            return d[3] + 0.1*(d[4]-d[3])*(t-30.);
        else if ( t<50. )
            return d[4] + 0.1*(d[5]-d[4])*(t-40.);
        else if ( t<60. )
            return d[5] + 0.1*(d[6]-d[5])*(t-50.);
        else
            return d[6];
    }
    
    // eq 9.27 de Olhoeft (1981)
    // c concentration en mole/l
    template<typename T>
    T corr_e_s(const T c) {
        return -13.*c + 1.065*c*c - 0.03006*c*c*c;
    }
    
    template<typename T>
    std::complex<T> ColeCole(const T es, const T ei, const T tr, const T f,
                             const T q, const T s) {
        // es: epsilon statique
        // ei: epsilon optique
        // tr: temps de relaxation (s)
        // f: frequence Hz
        // q = 1-alpha
        // s: conductivite
        T w = 2.*Utils_bg::pi*f;
        std::complex<T> tmp(ei, s/w);
		T tmp2 = 1.;
        return tmp + (es-ei)/(tmp2-pow(std::complex<T>(0.,w*tr), q));
    }
    
    template<typename T>
    std::complex<T> ColeCole_w(const T f, const T s, const T t, const T c) {
        // f -> frequence Hz
        // s -> conductivite S/m
        // t -> temperature C
		// c -> concentration en mole/l
        T e_s = (MalmbergMaryott(t) + corr_e_s(c)) * Utils_bg::epsilon0;
        T e_i = relax_w_epsilon_i(t)*Utils_bg::epsilon0;
        T tau_w = relax_w_tau(t);
        T q_w = 1.-relax_w_alpha(t);
        return ColeCole(e_s, e_i, tau_w, f, q_w, s);
    }
    
    template<typename T>
    std::complex<T> MaxwellWagner(const T f,const T delta,const T f1,const T f2){
        // eq  11 de druchinin00
        T tmp = delta/log(f2/f1);
        return std::complex<T>(tmp*(0.5*log((f*f+f2*f2)/(f*f+f1*f1))),
                               tmp*(atan(f2/f) - atan(f1/f)));
    }
    
    
    template<typename T>
    std::complex<T> permitArgile(const T e, const T t, const T wa, const T sw) {
        // e -> permittivite du grain d'argile
        // t -> temperature C
        // wa -> faction volumique d'eau liee au sein de l'argile
        // sw -> conductivite de l'eau liee
        
        // a 10 deg C, tableau 1 de Druchinin
        T a = 0.25;
        T f_bw = 1./relax_w_tau(10.0)/0.95;
        T f1_10 = f_bw/(1.+a);
        T f2_10 = f_bw*(1.+a);
        T de = 0.74*( MalmbergMaryott(10.0)-relax_w_epsilon_i(10.0) );
    }
    
    //
    // Classe DepolarizationFactor
    //
    template<typename T>
    class DepolarizationFactor{
    public:
        DepolarizationFactor(const std::valarray<T>);
        DepolarizationFactor(const T);
        
        T operator[](int n) { return (n>=0 && n<3 && valide) ? N[n] : 0.; }
        
    private:
        // 	std::valarray<T> N;
        T N[3];
        bool valide;
        
        void calculN(std::valarray<T>);
    };
    
    template<typename T>
        DepolarizationFactor<T>::DepolarizationFactor(const std::valarray<T> a) {
            if ( a.size() != 3 ) {
                valide = false;
                std::cerr << "Erreur: construction de DepolarizationFactor, argument invalide" << std::endl;
            }
            else {
                valide = true;
                // 	  N.resize(3);
            }
            this->calculN(a);
        }
    
    template<typename T>
        DepolarizationFactor<T>::DepolarizationFactor(const T xi) {
            valide = true;
            // 	N.resize(3);
            std::valarray<T> a(3);
            a[0] = 1.0;   a[1] = xi;   a[2] = xi;
            this->calculN(a);
        }
    
    
    template<typename T>
        void DepolarizationFactor<T>::calculN(const std::valarray<T> a) {
            // sihvola00, section 3.2.1
            if ( a[0]==a[1] && a[0]==a[2] ) {
                N[0] = N[1] = N[2] = 1.0/3.0;
            }
            else if (a[0]==a[1]) {
                if ( a[0] > a[2] ) {
                    double e = sqrt( (a[0]*a[0])/(a[2]*a[2]) - 1. );
                    N[2] = (1.+e*e)/(e*e*e)*(e-atan(e));
                    N[0] = N[1] = 0.5*(1.-N[2]);
                }
                else {
                    double e = sqrt( 1. - (a[1]*a[1])/(a[0]*a[0]) );
                    N[0] = (1-e*e)/(2.*e*e*e)*(log((1.+e)/(1-e) - 2.*e));
                    N[1] = N[2] = 0.5*(1.-N[0]);
                }
            }
            else if (a[0]==a[2]) {
                if ( a[0] > a[1] ) {
                    double e = sqrt( (a[0]*a[0])/(a[1]*a[1]) - 1. );
                    N[1] = (1.+e*e)/(e*e*e)*(e-atan(e));
                    N[0] = N[2] = 0.5*(1.-N[1]);
                }
                else {
                    double e = sqrt( 1. - (a[2]*a[2])/(a[0]*a[0]) );
                    N[0] = (1-e*e)/(2.*e*e*e)*(log((1.+e)/(1-e) - 2.*e));
                    N[2] = N[1] = 0.5*(1.-N[0]);
                }
            }
            else if (a[1]==a[2]) {
                if ( a[1] > a[0] ) {
                    double e = sqrt( (a[1]*a[1])/(a[0]*a[0]) - 1. );
                    N[0] = (1.+e*e)/(e*e*e)*(e-atan(e));
                    N[1] = N[2] = 0.5*(1.-N[0]);
                }
                else {
                    double e = sqrt( 1. - (a[2]*a[2])/(a[1]*a[1]) );
                    N[1] = (1-e*e)/(2.*e*e*e)*(log((1.+e)/(1-e) - 2.*e));
                    N[2] = N[0] = 0.5*(1.-N[1]);
                }
            }
        }
    
    
    //
    // Classe MaxwellGarnett
    //
    template<typename T>
    class MaxwellGarnett {
    public:
		// xi  -> ratio des axes de l'ellipsoide (Druchinin)
		MaxwellGarnett(const T xi=35.)
		: N(DepolarizationFactor<T>(xi)), c1(1.), c3(3.) {}
		MaxwellGarnett(const std::valarray<T>& xi)
		: N(DepolarizationFactor<T>(xi)), c1(1.), c3(3.) {}
		
		std::complex<T> operator()(const std::complex<T> e_i,
								   const std::complex<T> e_m, const T phi_m) {
			//
			// Maxwell-Garnett (formule 3.28 de Sihvola)
			//
			// e_i -> inclusion
			// e_m -> medium
			// e_e -> melange (e effectif)
			// phi_m -> fraction volumique du medium
			T f = 1.-phi_m;
			std::complex<T> e_e = e_m;
			std::complex<T> tmp;
			tmp = 0.;
			for (int n=0; n<3; ++n) tmp += (e_i-e_m)/(e_m + N[n]*(e_i-e_m));
			e_e *= f/c3 * tmp;
			tmp = 0.;
			for (int n=0; n<3; ++n) tmp += N[n]*(e_i-e_m)/(e_m + N[n]*(e_i-e_m));
			e_e /= c1 - f/c3 * tmp;
			return e_e+e_m;
		}
		
    private:
		typename DepolarizationFactor<T>::DepolarizationFactor N;
		const T c1;
		const T c3;
	};
    
    
    
    //
    // Classe HanaiBruggeman
    //
    template<typename T>
    class HanaiBruggeman {
    public:
        HanaiBruggeman() : p(std::valarray<T>(8)), xi(std::valarray<T>(4)), ftol(1.e-6)
	    {
            xi = 0.;  xi[0] = 1.;  xi[3] = 1.;
        }
        
        // e_i -> inclusion
        // e_m -> medium
        // phi_m -> fraction volumique du medium (egal a 1-phi_i)
        std::complex<T> operator()(const std::complex<T> e_i,
                                   const std::complex<T> e_m,
                                   const T W, const T phi_m)
	    {
            p[0] = 0.5*(e_i.real()+e_m.real());
            p[1] = 0.5*(e_i.imag()+e_m.imag());
            p[2] = e_i.real();	p[3] = e_i.imag();
            p[4] = e_m.real();	p[5] = e_m.imag();
            p[6] = W;           p[7] = phi_m;
 
//            clock_t t1 = clock();
//            if ( t1 == clock_t(-1) ) {
//                std::cerr << "pas de cloque...\n";
//                exit(1);
//            }
            Powell<T>::powell(p, xi, ftol, iter, fret, this->HanaiBruggemanDiff3);
//            clock_t t2 = clock();

//            p[0] = 0.5*(e_i.real()+e_m.real());
//            p[1] = 0.5*(e_i.imag()+e_m.imag());
//            clock_t t3 = clock();
//            Powell<T>::powell(p, xi, ftol, iter, fret, this->HanaiBruggemanDiff3);
//            clock_t t4 = clock();

            
//            std::cout << "STL, temps = " << 1000.0*double(t2-t1)/CLOCKS_PER_SEC << '\n';
//            std::cout << "BG, temps = " << 1000.0*double(t4-t3)/CLOCKS_PER_SEC << '\n';

            return std::complex<T>(p[0], p[1]);
        }
        
    private:
		std::valarray<T> p;
		std::valarray<T> xi;
		T ftol;
		long iter;
		T fret;
            
		static std::complex<T> HanaiBruggemanDiff(const std::complex<T> e_e,
												  const std::complex<T> e_i,
												  const std::complex<T> e_m,
												  const T W, const T phi_m)
		{
			// e_i -> inclusion
			// e_m -> medium
			// e_e -> melange (e effectif)
			// phi_m -> fraction volumique du medium (egal a 1-phi_i)
			return phi_m-((e_i-e_e)/(e_i-e_m))*pow(e_m/e_e, W);
		}

        static T HanaiBruggemanDiff2(const std::valarray<T>& x)
        {
			return abs( HanaiBruggemanDiff( std::complex<T>(x[0], x[1]),
                                            std::complex<T>(x[2], x[3]),
											std::complex<T>(x[4], x[5]),
											x[6], x[7]) );
		}
        
		static T HanaiBruggemanDiff3(const std::valarray<T>& x)
		{
            // e_e   -> x[0], x[1]
            // e_i   -> x[2], x[3]
            // e_m   -> x[4], x[5]
            // W     -> x[6]
            // phi_m -> x[7]
            
            // e_i - e_e -> a,b
            T a = x[2]-x[0];
            T b = x[3]-x[1];
            // e_i - e_m -> c,d
            T c = x[2]-x[4];
            T d = x[3]-x[5];
            
            // (e_i-e_e)/(e_i-e_m) -> e,f
            T tmp = c*c + d*d;
            T e = (a*c + b*d)/tmp;
            T f = (b*c - a*d)/tmp;
            
            // e_m/e_e -> a,b
            tmp = x[0]*x[0] + x[1]*x[1];
            a = (x[4]*x[0] + x[5]*x[1])/tmp;
            b = (x[5]*x[0] - x[4]*x[1])/tmp;
            
            // pow(e_m/e_e, W) -> c,d
            T arg = std::atan2(b, a);
            T g = std::pow(a*a + b*b, 0.5*x[6]);
            c = g * std::cos(x[6]*arg);
            d = g * std::sin(x[6]*arg);
            
            // ((e_i-e_e)/(e_i-e_m))*pow(e_m/e_e, W)  -> a,b
            a = e*c - f*d;
            b = e*d + f*c;
            
            a -= x[7];
            
            T mag = 1.;
            if ( x[0]<0. || x[1] < 0. ) mag = 1000.;
            
            return mag*sqrt(a*a + b*b);
        }
	};
    

    
    
    //
    // Classe HanaiBruggeman pour réels
    //
    template<typename T>
    class HanaiBruggemanR {
    public:
        HanaiBruggemanR() : p(std::valarray<T>(4)), ftol(1.e-6)
	    {
        }
        
        // e_i -> inclusion
        // e_m -> medium
        // phi_m -> fraction volumique du medium (egal a 1-phi_i)
        T operator()(const T e_i, const T e_m,
                     const T W, const T phi_m)
	    {
            T e_e;
            
            T bx = 0.5*(e_i+e_m);
            
            p[0] = e_i;
            p[1] = e_m;
            p[2] = W;
            p[3] = phi_m;
            
            Brent<T>::brent(e_i, bx, e_m, this->HanaiBruggemanDiff, ftol, e_e, p);
            return e_e;
        }
        
    private:
		std::valarray<T> p;
		T ftol;
        
		static T HanaiBruggemanDiff(const T e_e, const T e_i, const T e_m,
                                    const T W, const T phi_m)
		{
			// e_i -> inclusion
			// e_m -> medium
			// e_e -> melange (e effectif)
			// phi_m -> fraction volumique du medium (egal a 1-phi_i)
			return fabs( phi_m-((e_i-e_e)/(e_i-e_m))*pow(e_m/e_e, W) );
		}
        
	};
    
    //
    // Classe CRIM
    //
    template<typename T>
    class CRIM {
    public:
        CRIM() {}
        
        std::complex<T> operator()(const std::complex<T> e_i,
                                   const std::complex<T> e_m,
                                   const T phi_m) {
            std::complex<T> t = phi_m*sqrt(e_m) + (1.-phi_m)*sqrt(e_i);
            return t*t;
        }
    };
	
	//
	// Classe Chang & Liao
	//
	template<typename T>
	class ChangLiao {
    public:
		T operator() (const std::vector<T> &mu, const std::vector<T> &phi) {
			// eqn 23
			T somme = 0.0;
			for ( size_t n=0; n<mu.size(); ++n ) {
				somme += phi[n]*(mu[n]-1.)/(-2.*mu[n]+5.);
			}
			return (5.*somme + 1.)/(1. + 2.*somme);
		}
	};
    

	//
	// Classe Bussian solved by the method of Glover et al. 2010
	//
	template<typename T>
	class BussianGlover {
	public:
		BussianGlover() : p(std::valarray<T>(8)), xi(std::valarray<T>(4)), ftol(1.e-6)
		{
			xi = 0.;  xi[0] = 1.;  xi[3] = 1.;
		}
		
		// e_i -> inclusion
		// e_m -> medium
		// phi_m -> fraction volumique du medium (egal a 1-phi_i)
		std::complex<T> operator()(const std::complex<T> e_i,
								   const std::complex<T> e_m,
								   const T m, const T phi)
		{
			p[0] = 0.5*(e_i.real()+e_m.real());
			p[1] = 0.5*(e_i.imag()+e_m.imag());
			p[2] = e_i.real();	p[3] = e_i.imag();
			p[4] = e_m.real();	p[5] = e_m.imag();
			p[6] = m;           p[7] = phi;
			
			Powell<T>::powell(p, xi, ftol, iter, fret, this->BussianGloverDiff2);
			return std::complex<T>(p[0], p[1]);
		}
		
	private:
		std::valarray<T> p;
		std::valarray<T> xi;
		T ftol;
		long iter;
		T fret;
		
		static std::complex<T> BussianGloverDiff(const std::complex<T> e_e,
												 const std::complex<T> e_i,
												 const std::complex<T> e_m,
												 const T m, const T phi)
		{
			// e_i -> inclusion
			// e_m -> medium
			// e_e -> melange (e effectif)
			// m  cementation
			// phi -> porosite
			return (e_e-1.)/pow(e_e, 1./m) - phi * (e_m/e_i-1.)/pow(e_m/e_i, 1./m);
		}
		
		static T BussianGloverDiff2(const std::valarray<T>& x)
		{
			return abs( BussianGloverDiff(std::complex<T>(x[0], x[1]),
										  std::complex<T>(x[2], x[3]),
										  std::complex<T>(x[4], x[5]),
										  x[6], x[7]) );
		}
		
	};
	
};  // end namespace



#endif
