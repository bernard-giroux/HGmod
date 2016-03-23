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
#ifndef __PRIDE_H__
#define __PRIDE_H__

/*
 *
 *  Modèle de conductivité de Steve Pride
 *
 *
 *
 * Références
 *
 * @article{pride94,
 *   author = {Steve Pride},
 *   journal = {Physical Review B},
 *   number = {21},
 *   pages = {15678--15696},
 *   title = {Governing equations for the coupled electromagnetics and acoustic of porous media},
 *   volume = {50},
 *   year = {1994}
 * }
 *
 * @article{carcione03,
 *   author = {Jos\'e M. Carcione and G\'eza Seriani and Davide Gei},
 *   journal = {Journal of Applied Geophysics},
 *   pages = {177--191},
 *   title = {Acoustic and electromagnetic properties of soils saturated with salt water and {NAPL}},
 *   volume = {52},
 *   year = {2003}
 * }
 *
 */


#include <complex>
#include <cmath>

#include "constantes.h"
#include "Fluide.h"


template<typename T>
class Pride {
public:
	Pride(const Fluide<T> &fl, const T t=2.5, const T w=1.0)
    : fluide(fl), iT(1./t), omega(w)
    {
        zeta = calculZeta( fluide.getSalinite() );
        update_P();
        c = Cem_ReCos(w);
    }
	
	T sigma(const T porosite, const T Lambda) const
	{
		return porosite*iT*(fluide.getSigma() + 2.*c/Lambda);
	}
	
	T sigma(const T porosite, const T Lambda, const T Saturation, const T n) const
	{
        if ( Saturation < 0.999 )
            return porosite*iT*(fluide.getSigma()*std::pow(Saturation,n) + 2.*c/Lambda);
        else
            return porosite*iT*(fluide.getSigma() + 2.*c/Lambda);
	}
	
	T sigma(const T porosite, const T Lambda, const T Saturation, const T n, const Fluide<T>& fl)
	{
        zeta = calculZeta( fl.getSalinite() );
        update_P();
        if ( Saturation < 0.999 )
            return porosite*iT*(fl.getSigma()*std::pow(Saturation,n) +
                                2.*Cem_ReCos(omega)/Lambda);
        else
            return porosite*iT*(fl.getSigma() +
                                2.*Cem_ReCos(omega)/Lambda);
	}
	
	T sigmaHF(T porosite, T Lambda, T omega)
	{
		T c2 = Cem_ReCos(omega);
		return porosite*iT*(fluide.getSigma() + 2.*c2/Lambda);
	}
	
	T sigmaHF(T porosite, T Lambda, T omega, T Saturation, const T n)
	{
		T c2 = Cem_ReCos(omega);
        if ( Saturation < 0.999 )
            return porosite*iT*(fluide.getSigma()*std::pow(Saturation,n) + 2.*c2/Lambda);
        else
            return porosite*iT*(fluide.getSigma() + 2.*c2/Lambda);
	}
	
	T getTortuosite() const { return 1./iT; }
	
	static T ppt2mol_l(T ppt, T rho=1.) { return ppt*rho/58.443; }
	static T ppm2mol_l(T ppm, T rho=1.) { return ppm*rho/58443.; }
	
	// conversion Surface spécifique [m2/m3] à Lambda [m3/m2]
	static T S2Lambda(T S, T phi) {	return 2.*phi/(S*(1.-phi)); }
	
	// conversion perméabilité [mD] à Lambda [m3/m2]
	static T k2Lambda(T k, T phi, T t=2.5, T E=10.)
	{
		return sqrt(E*t*(k*9.8692327e-16)/phi);
	}
	
	// conversion Surface spécifique [m2/m3] à perméabilité [mD]
	// d'après l'eq B-19 de Carcione 
	static T S2k(T S, T phi, T t=2.5, T E=10.)
	{
		T L =  2.*phi/(S*(1.-phi));
		return 1.013249996628410e+15*L*L*phi/(E*t);
	}
	
    // equation 54, Pride et Morgan, 1991
    static T calculZeta(T s) { return 0.008 + 0.026*std::log10(s); }

    
private:
	const Fluide<T> fluide;
	const T iT;          // 1/tortuosité
	const T omega;       // frequence angulaire
	T c;
    T P;
    T zeta;              // potentiel zeta [V]
    
    // eq 194 de Pride
    T Cem() const
    {
        const T* ez = fluide.getCharge();
        const T* N = fluide.getConcentration();
        const T* b = fluide.getMobilite();
        T d = fluide.getEpaisseurDebye();
        
        T tmp = ez[0]*ez[0]*b[0]*N[0]*
        (std::exp((-ez[0]*zeta)/(2*Utils_bg::k*fluide.getTemperature()))-1.0);
        tmp += ez[1]*ez[1]*b[1]*N[1]*
        (std::exp((-ez[1]*zeta)/(2*Utils_bg::k*fluide.getTemperature()))-1.0);
        return 2.0*d*tmp;
    }
    
    // eq 206 de Pride
    std::complex<T> Cos(T w) const
    {
        // à la freq angulaire w
        T tmp = P*(fluide.getEpsilon()*fluide.getEpsilon()*zeta*zeta)/
        (2.0*fluide.getEpaisseurDebye()*fluide.getViscosite());
        return tmp/(1.0-(2.0*Utils_bg::i32*fluide.getEpaisseurDebye())/(P*delta(w)));
    }
    
    T Cem_ReCos( T w ) const
    {
        T tmp1 = this->Cem();
        std::complex<T> tmp2 = this->Cos(w);
        return tmp1+tmp2.real();
    }

    // eq 207 de Pride
    void update_P()
    {
        const T* ez = fluide.getCharge();
        const T* N = fluide.getConcentration();
        T d = fluide.getEpaisseurDebye();

        P  = N[0]*(std::exp((-ez[0]*zeta)/(2.0*Utils_bg::k*fluide.getTemperature()))-1.0);
        P += N[1]*(std::exp((-ez[1]*zeta)/(2.0*Utils_bg::k*fluide.getTemperature()))-1.0);
        P *= (8.*Utils_bg::k*fluide.getTemperature()*d*d)/(fluide.getEpsilon()*zeta*zeta);
    }
    
    
    // delta: 0.01 pour passer de cp à p, et 0.1 pour passer de g/cm à kg/m
    T delta(T w) const 
    {
        return std::sqrt(0.001*fluide.getViscosite()/(w*fluide.getDensite()));
    }

};


template<typename T>
class WaxmanSmits {
public:
	WaxmanSmits(T ae, T s, T nn=2., T mm=1.37, T aa=0.88)
	: alpha_e(ae), sigma_f(s), n(nn), m(mm), a(aa)
	{
		calculB();
	}

	T sigma(T porosite, T Surf_specifique, T Saturation)
	{
		calculQv( porosite, Surf_specifique );
		if ( Saturation > 0.999 ) 
			return a/std::pow(porosite, m) * ( sigma_f + B*Qv );
		else
			return std::pow(Saturation, n)*a/std::pow(porosite, m) *
				( sigma_f + B*Qv/Saturation );
	}

private:
  T alpha_e;    // densite de charge surfacique [meq/m2]
  T sigma_f;    // conductivité du fluide
  T n;
  T m;          // facteur de cimentation
  T a;
  T B;          // counterion mobility
  T Qv;         // shalyness parameter

  void calculB() { B = 3.83*(1.-0.83*std::exp(-0.5*sigma_f)); } // eq 9-111 dans Schon
  void calculQv(T phi, T Surf_specifique)
  { // Surf_specifique en m2/m3
	  Qv = 1.e-6*alpha_e*Surf_specifique*(1.-phi)/phi;
  }
  
};

template<typename T>
class WaxmanSmits2 {
public:
	WaxmanSmits2(T d, T ae, T s, T nn=2., T mm=1.37, T aa=0.88)
	: densite(d), alpha_e(ae), sigma_f(s), n(nn), m(mm), a(aa)
	{
		calculB();
	}
	
	T sigma(T porosite, T Surf_specifique, T Saturation)
	{
		calculCEC( Surf_specifique );
		calculQv( porosite );
		if ( Saturation > 0.9999 ) 
			return a/std::pow(porosite, m) * ( sigma_f + B*Qv );
		else
			return std::pow(Saturation, n)*a/std::pow(porosite, m) *
				( sigma_f + B*Qv/Saturation );
	}

private:
	T densite;    // densite dans grains solides [g/cm3]
	T alpha_e;    // densite de charge surfacique [meq/m2]
	T sigma_f;    // conductivité du fluide
	T n;
	T m;          // facteur de cimentation
	T a;
	T B;          // counterion mobility
	T Qv;         // shalyness parameter
	T CEC;        // cation exchange capacity [meq/100g]
	
	void calculB() { B = 3.83*(1.-0.83*std::exp(-0.5*sigma_f)); }
	void calculCEC(T Surf_specifique) { // Surf_specifique en m2/100g
		CEC = alpha_e*Surf_specifique; }
	void calculQv(T phi) { Qv = 0.01*densite*CEC*(1.-phi)/phi; }
	
};

#endif
