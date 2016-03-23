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
 * }
 * 
 * @BOOK{weast80,
 *	 title = {Handbook of Chemistry and Physics},
 *	 publisher = {CRC Press},
 *	 year = {1980},
 *	 author = {R.C. Weast and M. J. Astle},
 *	 address = {Boca Raton, Florida},
 *	 edition = {61},
 *	 owner = {giroux},
 * }
 *
 */

#ifndef __MILIEUPOREUX_H__
#define __MILIEUPOREUX_H__

#include <ostream>

#include "constantes.h"

#include "structHGmod.h"
#include "Covariance.h"
#include "Generateur.h"
#include "Granulo.h"
#include "Lima.h"
#include "ModelesMelange.h"
#include "Pride.h"
#include "Raymer.h"
#include "Krigeage.h"
#include "Messages.h"

extern Messages msg;

using namespace ModelesMelange;

template<typename T, template<typename>class GEN, typename COV>
class MilieuPoreux
{
public:
	MilieuPoreux(const paramPhysique<T>&, const Generateur<COV>*);
	//  MilieuPoreux(const MilieuPoreux&);
	
	~MilieuPoreux()
	{
		delete fluide; delete pride; delete lima;
        delete granulo; delete sCov; delete sKrige;
	}
	    
    bool conditionne() const
    {
        return pPhys.pGen.pCond.conditionne;
    }
    
    const paramPhysique<T>& get_paramPhysique() const { return pPhys; }
    
    double get_covariance(const std::valarray<double> x1,
                          const std::valarray<double> x2) const
    {
        return generateur->get_covariance(x1, x2);
    }
    
    double getPorositeMoyenne() const
    {
        return pPhys.mPorosite;
    }
	
//	T getSigma(const double x, const double y, const double z) const
//	{
//		return calculSigma(x, y, z);
//	}
//	
//	T getSigma(const double x, const double y, const double z, const double phi) const
//	{
//		return calculSigma(x, y, z, phi);
//	}
	
	T getSigma(const double S, const double phi) const
	{
		return calculSigma(S, phi);
	}
    
	T getSigma(const double S, const double phi, const Fluide<T>& fl) const
	{
		return calculSigma(S, phi, fl);
	}

    std::complex<T> getSigmaCplx(const double S, const double phi, const double f) const
	{
		return calculSigmaCplx(S, phi, f);
	}
    
	std::complex<T> getSigmaCplx(const double S, const double phi, const double f,
                                 const Fluide<T>& fl) const
	{
		return calculSigmaCplx(S, phi, f, fl);
	}
    
    T getChargeabilite(const double S, const double phi) const
    {
        return 0.0;
    }        
	
//	std::complex<T> getPermit(const double x, const double y,
//							  const double z, const double f) const
//	{
//		return calculPermit(x, y, z, f);
//	}
//	
//	std::complex<T> getPermit(const double x, const double y,
//							  const double z, const double f,
//                              const double phi) const
//	{
//		return calculPermit(x, y, z, f, phi);
//	}
	
	std::complex<T> getPermit(const double S, const double phi, const double f) const
	{
		return calculPermit(S, f, phi);
	}
	
	std::complex<T> getPermit(const double S, const double phi, const double f,
                              const Fluide<T>& fl) const
	{
		return calculPermit(S, f, phi, fl);
	}
	
	T getPermea(const double x, const double y, const double z) const
	{
		return pPhys.permeabilite;
	}
	
	T getSigmaM(const double x, const double y, const double z) const
	{
		return pPhys.stochastique ? calculSigmaM(x, y, z) : pPhys.conductivite;
	}
	
	std::complex<T> getPermitM(const double x, const double y,
							   const double z) const
	{
		return pPhys.stochastique ? calculPermitM(x, y, z, pPhys.f) :
		std::complex<T>(pPhys.permittivite, 0.0);
	}
	
	T getPermeaM(const double x, const double y, const double z) const
	{
		return pPhys.permeabilite;
	}
	
	T getPorositeTotale(const double x,const double y,const double z) const
	{
		return pPhys.mPorosite + generateur->z(x, y, z);
	}
	
	T getSaturation(const double x, const double y, const double z) const
	{
		return calculSaturation( pPhys.pSat, z );
	}
	
	T getFractionArgile(const double x, const double y, const double z)	const
	{
		return getFractionArgile( getPorositeTotale(x,y,z) );
	}
	
	T getFractionArgile(const double phi) const
	{
        if (pPhys.phi_c<0.0) return 0.0;
		T f = (phi-pPhys.phi_s)/(pPhys.phi_c-1.);
        if ( phi > pPhys.phi_s ) f = 0.;
        if ( phi < pPhys.phi_s*pPhys.phi_c ) f = pPhys.phi_s;
		return f;
	}
	
	// retourne la perméabilité hydraulique en mD
	T getPermeaHydro(const double x, const double y, const double z) const
	{
		return getPermeaHydro(getPorositeTotale(x,y,z));
	}
	
	// retourne la perméabilité hydraulique en mD
	T getPermeaHydro(const double phi) const
	{
		T c = getFractionArgile(phi);
		c *= (1.-pPhys.phi_c) / (1.- phi);  // correction term published in errata
		T Ss = pPhys.S_c*c + (1.-c)*pPhys.S_s;
        return pPhys.KozenyCarman ? KozenyCarman(phi, Ss, pPhys.T, pPhys.facteurPermea) :
            Pride<T>::S2k(Ss, phi, pPhys.T, pPhys.facteurPermea);  // modèle de Pride
	}
	
	T getFreqTravail() const { return pPhys.f; }
	
	double getVariableAleatoire(const double x,const double y,const double z)
		const { return generateur->z(x, y, z); }
	
	void affiche_param(std::ostream&);
	
	T getSalinite(const std::valarray<double>& x) const {
        return sKrige==0 ? pPhys.TDS : sKrige->Z_est(x);
    }
	T getTemperature() const { return pPhys.temperature; }
    
	const GEN<COV>* getGenerateur() const { return generateur; }
    
    T getVp(const T phi) const { return raymer->V(phi); }
private:
	const paramPhysique<T> pPhys;
	const GEN<COV>* generateur;
	const T iepsilon0;
	
	Fluide<T>* fluide;
	Pride<T>* pride;
    Lima<T>* lima;
    Raymer<T>* raymer;
	Granulo* granulo;
	COV* sCov;
	KrigeageOrdinaire<COV>* sKrige;

	T sigma_c;   // conductivité de l'argile (S/m)
    T sw;        // conductivité de l'eau (S/m)
		
	T calculPorosite(const double, const double, const double) const;
//	T calculSigma(const double, const double, const double) const;
//	T calculSigma(const double, const double, const double, const double) const;
	T calculSigma(const double, const double) const;
	T calculSigma(const double, const double, const Fluide<T>&) const;
    std::complex<T> calculSigmaCplx(const double, const double, const double) const;
    std::complex<T> calculSigmaCplx(const double, const double, const double,
                                    const Fluide<T>&) const;
//	std::complex<T> calculPermit(const double, const double, const double,
//								 const T) const;
//	std::complex<T> calculPermit(const double, const double, const double,
//								 const T, const double) const;
	std::complex<T> calculPermit(const double, const T, const double) const;
	std::complex<T> calculPermit(const double, const T, const double, const Fluide<T>&) const;
	T calculSigmaM(const double, const double, const double) const;
	T calculPermitM(const double, const double, const double) const;
	std::complex<T> calculPermitM(const double, const double, const double,
								  const T) const;
	T calculSaturation( const paramSaturation<T>& p, const double z ) const;
	T calculSigmaEau(const T) const;
	
	std::complex<T> eq_3(const T, const T, const T, const T) const;
	std::complex<T> eq_3(const T, const T, const T, const T, const Fluide<T>&) const;
	std::complex<T> eq_4(const T, const T, const T, const T) const;
	std::complex<T> eq_4(const T, const T, const T, const T, const Fluide<T>&) const;
	std::complex<T> eq_5(const T, const T, const T, const T, const T) const;
	std::complex<T> eq_5(const T, const T, const T, const T, const T, const Fluide<T>&) const;
	std::complex<T> eq_6(const T, const T, const T, const T, const T) const;
	std::complex<T> eq_6(const T, const T, const T, const T, const T, const Fluide<T>&) const;
	
	std::complex<T> druchinin(const T phi_c, const T phi_s,
							  const T phi_w, const T phi_a,
							  const T f) const
	{
		if ( pPhys.S_c == pPhys.S_s ) {
			if ( phi_w < 0.11 ) return eq_3(phi_c+phi_s, phi_w, phi_a, f);
			else return eq_4(phi_c+phi_s, phi_w, phi_a, f);
		}
		if ( phi_c+phi_s < phi_w ) return eq_5(phi_c, phi_s, phi_w, phi_a, f);
		else return eq_6(phi_c, phi_s, phi_w, phi_a, f);
	}
	
	std::complex<T> druchinin(const T phi_c, const T phi_s,
							  const T phi_w, const T phi_a,
							  const T f, const Fluide<T>& fl) const
	{
		if ( pPhys.S_c == pPhys.S_s ) {
			if ( phi_w < 0.11 ) return eq_3(phi_c+phi_s, phi_w, phi_a, f, fl);
			else return eq_4(phi_c+phi_s, phi_w, phi_a, f, fl);			
		}
		if ( phi_c+phi_s < phi_w ) return eq_5(phi_c, phi_s, phi_w, phi_a, f, fl);
		else return eq_6(phi_c, phi_s, phi_w, phi_a, f, fl);
	}
	
};



template<typename T, template<typename>class GEN, typename COV>
MilieuPoreux<T,GEN,COV>::MilieuPoreux(const paramPhysique<T>& par,
									  const Generateur<COV>* g)
: pPhys(par), generateur(g), iepsilon0(1./Utils_bg::epsilon0), fluide(0),
	pride(0), lima(0), raymer(0), granulo(0), sCov(0), sKrige(0)
{
	T s = Pride<T>::ppm2mol_l(pPhys.TDS);
	T z[] = {1., -1.};
	T R[] = {116.e-12, 167.e-12}; // rayon ionique en m pour Na et Cl respectivement
	T eta = Utils_bg::viscositeEau(pPhys.temperature);  // viscosité en cP
	
    // mobilité des ions en (m/s)/N (facteur 1000 pour passer de cP à kg/m.s)
	T b[] = {1000./(6.*Utils_bg::pi*eta*R[0]), 1000./(6.*Utils_bg::pi*eta*R[1])};
    
	T rho_w = Utils_bg::masseVolEau(pPhys.TDS, pPhys.temperature);  // densité de l'eau
	T epsilon = (MalmbergMaryott( pPhys.temperature ) + corr_e_s(s))*Utils_bg::epsilon0;
    T pH = 7.0;
	fluide = new Fluide<T>(s, 273.+pPhys.temperature, rho_w, eta, epsilon, pH, z, b);
	
	sw = fluide->getSigma();  // conductivité de l'eau libre (S/m)
	T omega=2.*Utils_bg::pi;  // 1 Hz
	pride = new Pride<T>( *fluide, pPhys.T, omega );
    
	granulo = new Granulo(pPhys.granuloDiam, pPhys.granuloPassant);
    if ( pPhys.granuloMasseVol.size() > 0 ) {
        granulo->setMasseVolumique( pPhys.granuloMasseVol );
    }
    if ( pPhys.granuloFactDepol.size() > 0 ) {
        granulo->setFacteurDepolarisation( pPhys.granuloFactDepol );
    }
    
    lima = new Lima<T>( fluide );
    
    raymer = new Raymer<T>( pPhys.V0, pPhys.Vfl );
    
	if ( pPhys.sigma_c < 0. )  // si on a donné une conductivité négative pour l'argile
    {
		T Lambda = Pride<T>::S2Lambda(pPhys.S_c, pPhys.phi_c);
		sigma_c = pride->sigma(pPhys.phi_c, Lambda);
	}
	else
		sigma_c = pPhys.sigma_c;
	
	if ( pPhys.sData.conditionne ) {
		// la salinité n'est pas constante
		sCov = new COV();
		sCov->remplire(pPhys.sCov);
		sKrige = new KrigeageOrdinaire<COV>(sCov, &(pPhys.sData));
	}
	else {
		sCov = 0;
		sKrige = 0;
	}
}
    
template<typename T, template<typename>class GEN, typename COV>
inline T MilieuPoreux<T,GEN,COV>::calculPorosite(const double x, const double y,
												 const double z) const
{
    T phi = pPhys.mPorosite + generateur->z(x, y, z);
    return phi<0. ? 0. : (phi>1. ? 1. : phi);
}

template<typename T, template<typename>class GEN, typename COV>
T MilieuPoreux<T,GEN,COV>::calculSigma(const double S, const double phi) const
{
	T c = getFractionArgile(phi);
	c *= (1.-pPhys.phi_c) / (1.-phi);  // correction term published in errata
	T Ss = pPhys.S_c>0.0 ? pPhys.S_c*c + (1.-c)*pPhys.S_s : pPhys.S_s;
	T Lambda = Pride<T>::S2Lambda(Ss, phi);
	return pride->sigma(phi, Lambda, S, pPhys.n);
}

template<typename T, template<typename>class GEN, typename COV>
T MilieuPoreux<T,GEN,COV>::calculSigma(const double S, const double phi,
                                       const Fluide<T>& fl) const
{
	T c = getFractionArgile(phi);
	c *= (1.-pPhys.phi_c) / (1.-phi);  // correction term published in errata
	T Ss = pPhys.S_c>0.0 ? pPhys.S_c*c + (1.-c)*pPhys.S_s : pPhys.S_s;
	T Lambda = Pride<T>::S2Lambda(Ss, phi);
	return pride->sigma(phi, Lambda, S, pPhys.n, fl);
}

template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::calculSigmaCplx(const double S,
                                                         const double phi,
                                                         const double f) const
{
    std::vector<double> p = granulo->getFractionVolumique();  // fraction du volume des grains solides
    std::vector<double> W = granulo->getFacteurDepolarisation();
    for (size_t n=0; n<p.size(); ++n) {
        // pour obtenir la fraction du volume total
        p[n] *= (1.0-phi);
    }
    double omega = 2.*Utils_bg::pi*f;
    std::complex<T> s = fluide->getSigma(omega);
    T pm = phi*S;  // initialement, le "medium" est l'eau
    
    std::complex<T> s2;
    HanaiBruggeman<T> hb;
    for (ptrdiff_t n=p.size()-1; n>=0; --n) {  // on commence par les grains fins
        s2 = lima->sigma_fixman( granulo->getDiametre(n), omega);
        s = hb(s2, s, W[n], pm/(pm+p[n]));
        pm += p[n];  // la fraction volumique du milieu contien maintenant la dernière fraction solide
    }
    if ( S < 0.99) {
        // on tient compte de l'air
        pm = 1.0 - phi*(1.0-S);
        // 1/(1-W) = 1.5 pour W=1/3
        s *= std::pow(pm, 1.5);
    }
    return s;
}

template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::calculSigmaCplx(const double S,
                                                         const double phi,
                                                         const double f,
                                                         const Fluide<T>& fl) const
{
    std::vector<double> p = granulo->getFractionVolumique();  // fraction du volume des grains solides
    std::vector<double> W = granulo->getFacteurDepolarisation();
    for (size_t n=0; n<p.size(); ++n) {
        // pour obtenir la fraction du volume total
        p[n] *= (1.0-phi);
    }
    double omega = 2.*Utils_bg::pi*f;
    std::complex<T> s = fl.getSigma(omega);
    T pm = phi*S;  // initialement, le "medium" est l'eau
    
    std::complex<T> s2;
    HanaiBruggeman<T> hb;
    for (ptrdiff_t n=p.size()-1; n>=0; --n) {  // on commence par les grains fins
        s2 = lima->sigma_fixman( granulo->getDiametre(n), omega);
        s = hb(s2, s, W[n], pm/(pm+p[n]));
        pm += p[n];  // la fraction volumique du milieu contien maintenant la dernière fraction solide
    }
    if ( S < 0.99) {
        // on tient compte de l'air
        pm = 1.0 - phi*(1.0-S);
        // 1/(1-W) = 1.5 pour W=1/3
        s *= std::pow(pm, 1.5);
    }
    return s;
}

//template<typename T, template<typename>class GEN, typename COV>
//std::complex<T> MilieuPoreux<T,GEN,COV>::calculPermit(const double x,
//													  const double y,
//													  const double z,
//													  const T f) const
//{
//	return calculPermit(x, y, z, f, calculPorosite(x,y,z));
//}
//
//template<typename T, template<typename>class GEN, typename COV>
//std::complex<T> MilieuPoreux<T,GEN,COV>::calculPermit(const double x,
//													  const double y,
//													  const double z,
//													  const T f,
//													  const double phi) const
//{
//	T c = getFractionArgile(phi);
//	T phi_c = c;
//	T phi_w = phi*(1.-c)*getSaturation(x,y,z);
//	T phi_a = phi*(1.-c)*(1.-getSaturation(x,y,z));
//	T phi_s = 1.-phi_c-phi_w-phi_a;
//	
//	return druchinin(phi_c, phi_s, phi_w, phi_a, f);
//}

template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::calculPermit(const double S,
													  const T f,
													  const double phi) const
{
	T c = getFractionArgile(phi);
	T phi_c = c;
	T phi_w = phi*(1.-c)*S;
	T phi_a = phi*(1.-c)*(1.-S);
	T phi_s = 1.-phi_c-phi_w-phi_a;
	
	return druchinin(phi_c, phi_s, phi_w, phi_a, f);
}

template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::calculPermit(const double S,
													  const T f,
													  const double phi,
													  const Fluide<T>& fl) const
{
	T c = getFractionArgile(phi);
	T phi_c = c;
	T phi_w = phi*(1.-c)*S;
	T phi_a = phi*(1.-c)*(1.-S);
	T phi_s = 1.-phi_c-phi_w-phi_a;
	
	return druchinin(phi_c, phi_s, phi_w, phi_a, f, fl);
}

template<typename T, template<typename>class GEN, typename COV>
T MilieuPoreux<T,GEN,COV>::calculSigmaM(const double x, const double y,
										const double z) const
{
	return calculSigma(x, y, z, pPhys.mPorosite);     // porosite totale moyenne
}


template<typename T, template<typename>class GEN, typename COV>
T MilieuPoreux<T,GEN,COV>::calculPermitM(const double x, const double y,
										 const double z) const
{
	return calculPermit(x, y, z, fluide, pPhys.mPorosite);   // porosite totale moyenne
}

template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::calculPermitM(const double x,
													   const double y,
													   const double z,
													   const T f) const
{
	return calculPermit(x, y, z, f, pPhys.mPorosite);
}


template<typename T, template<typename>class GEN, typename COV>
T MilieuPoreux<T,GEN,COV>::calculSaturation( const paramSaturation<T>& p,
											 const double z ) const
{
	if ( z >= p.profNappe ) return p.Smax;
	if ( p.vG ) {  // van Genuchten 1980
		T h = p.profNappe - z;  // pressure head
		T m = p.vG_m;
		if ( p.vG_m == 0.0 ) m = 1.0 - 1.0/p.vG_n;
		return p.Smin + (p.Smax-p.Smin) * 1.0/pow( 1 + pow(p.vG_alpha*h, p.vG_n), m);
	}
	else {
		if ( z <= p.profNappe-p.epaisseurTransition ) return p.Smin;
		else {
			T x1 = (p.profNappe-z)/p.epaisseurTransition;
			T x = 2.5 - 5.*x1;
			//		T tmp = 1.75;
			//		return p.Smin + (1.-pow(x1,tmp))*(p.Smax-p.Smin)*
			//            (0.5+0.5*Utils_bg::bg_erf( x ) );
			return p.Smin + (p.Smax-p.Smin)*(0.5+0.5*Utils_bg::bg_erf( x ) );
		}
	}
}

template<typename T, template<typename>class GEN, typename COV>
inline T MilieuPoreux<T,GEN,COV>::calculSigmaEau(const T t) const
{
	// t -> temperature en C
	return 0.00016 * pPhys.TDS * ( 1.+ 0.02*(t-25.0) );
}

template<typename T, template<typename>class GEN, typename COV>
void MilieuPoreux<T,GEN,COV>::affiche_param(std::ostream& out)
{
	out << " - " << msg.getString("Propriétés") << '\n';
	out << "   - " << msg.getString("Porosité moyenne") << ": " << pPhys.mPorosite << '\n';
	out << "   - " << msg.getString("Température [°C]") << ": " << pPhys.temperature << '\n';
	out << "   - " << msg.getString("Permittivité relative du sable") << ": " << pPhys.Km << '\n';
	out << "   - " << msg.getString("Porosité du sable pur") << ": " << pPhys.phi_s << '\n';
	out << "   - " << msg.getString("Porosité de l'argile pure") << ": " << pPhys.phi_c << '\n';
	out << "   - " << msg.getString("Surface spécifique du sable pur [m²/m³]") << ": " << pPhys.S_s << '\n';
	out << "   - " << msg.getString("Surface spécifique de l'argile pure [m²/m³]") << ": " << pPhys.S_c << '\n';
	out << "   - " << msg.getString("Permittivité relative statique de l'argile") << ": " << pPhys.Ks_c << '\n';
	out << "   - " << msg.getString("Permittivité relative optique de l'argile") << ": " << pPhys.Ki_c << '\n';
	out << "   - " << msg.getString("Fréquence de relaxation de l'argile [Hz]") << ": " << pPhys.fr_c << '\n';
	out << "   - " << msg.getString("Exposant Cole-Cole - argile") << ": " << pPhys.ccq_c << '\n';
	out << "   - " << msg.getString("Tortuosité") << ": " << pPhys.T << '\n';
	if (pPhys.sData.conditionne) {
		out << "   - " << msg.getString("Modèle de salinité hétérogène") << '\n';
	} else {
		out << "   - " << msg.getString("Salinité [ppm] (homogène)") << ": " << pPhys.TDS << '\n';
	}
	if ( pPhys.pSat.vG ) {
		out << "   - " << msg.getString("Modèle de saturation de van Genuchten appliqué") << '\n';
		T m = pPhys.pSat.vG_m;
		if ( pPhys.pSat.vG_m == 0.0 ) m = 1.0 - 1.0/pPhys.pSat.vG_n;
		out << "     - " << "alpha: " << pPhys.pSat.vG_alpha << '\n';
		out << "     - " << "m: " << m << '\n';
		out << "     - " << "n: " << pPhys.pSat.vG_n << '\n';
	}
	if (generateur != 0)
        generateur->affiche_param( out );
		
}

//
// equation 3b de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_3(const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f // -> frequence
											  ) const
{
	std::complex<T> e_w;
	if ( f > 0.1 ) {
		e_w = ColeCole_w(f, sw, pPhys.temperature,
                         Pride<T>::ppm2mol_l(pPhys.TDS))*iepsilon0;
	}
	else {
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(Pride<T>::ppm2mol_l(pPhys.TDS));
	}
	// medium: sable, inclusion: eau
	T phi = phi_s+phi_w;
	MaxwellGarnett<T> mg(1.0);
	std::complex<T> e = mg( e_w, std::complex<T>(pPhys.Km, 0.0), phi_s/phi );
	// medium: e, inclusion: air
	if ( phi_a > 0.0 ) {
		HanaiBruggeman<T> hb;
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	}
	else
		return e;
}	
//
// equation 3b de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_3(const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f, // -> frequence
											  const Fluide<T>& fl) const
{
	std::complex<T> e_w;
	if ( f > 0.1 ) {
		e_w = ColeCole_w(f, fl.getSigma(), pPhys.temperature, fl.getSalinite())*iepsilon0;
	}
	else {
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(fl.getSalinite());
	}
	// medium: sable, inclusion: eau
	T phi = phi_s+phi_w;
	MaxwellGarnett<T> mg(1.0);
	std::complex<T> e = mg( e_w, std::complex<T>(pPhys.Km, 0.0), phi_s/phi );
	// medium: e, inclusion: air
	if ( phi_a > 0.0 ) {
		HanaiBruggeman<T> hb;
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	}
	else
		return e;
}	

//
// equation 4 de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_4(const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f // -> frequence
											  ) const
{
	std::complex<T> e_w;
	if ( f > 0.1 ) {
		e_w = ColeCole_w(f, sw, pPhys.temperature,
                         Pride<T>::ppm2mol_l(pPhys.TDS))*iepsilon0;
	}
	else {
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(Pride<T>::ppm2mol_l(pPhys.TDS));
	}
	// medium: sable, inclusion: eau
	T phi = phi_s+phi_w;
	HanaiBruggeman<T> hb;
	std::complex<T> e = hb( e_w, std::complex<T>(pPhys.Km, 0.0), 1./3., phi_s/phi );
	// medium: e, inclusion: air
	if ( phi_a > 0.0 )
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	else
		return e;
}	
//
// equation 4 de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_4(const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f, // -> frequence
											  const Fluide<T>& fl) const
{
	std::complex<T> e_w;
	if ( f > 0.1 ) {
		e_w = ColeCole_w(f, fl.getSigma(), pPhys.temperature, fl.getSalinite())*iepsilon0;	}
	else {
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(fl.getSalinite());
	}
	// medium: sable, inclusion: eau
	T phi = phi_s+phi_w;
	HanaiBruggeman<T> hb;
	std::complex<T> e = hb( e_w, std::complex<T>(pPhys.Km, 0.0), 1./3., phi_s/phi );
	// medium: e, inclusion: air
	if ( phi_a > 0.0 )
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	else
		return e;
}	

//
// equation 5 de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_5(const T phi_c, // -> argile
											  const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f // -> frequence
											  ) const
{
	std::complex<T> e_w, e_c;
	if ( f > 0.1 ) {
		e_w = ColeCole_w(f, sw, pPhys.temperature,
                         Pride<T>::ppm2mol_l(pPhys.TDS))*iepsilon0;
		T tmp = 1./pPhys.fr_c;
		e_c = ColeCole(pPhys.Ks_c, pPhys.Ki_c, tmp, f,
					   pPhys.ccq_c, sigma_c*iepsilon0);
	}
	else {
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(Pride<T>::ppm2mol_l(pPhys.TDS));
		e_c = pPhys.Ks_c;
	}
	
	// medium: eau, inclusion: argile
	T phi = phi_w+phi_c;
	MaxwellGarnett<T> mg(35.0);
	std::complex<T> e = mg( e_c, e_w, phi_w/phi );
	// medium: sable, inclusion: e
	phi += phi_s;
	e = mg(e, std::complex<T>(pPhys.Km, 0.0), phi_s/phi);
	// medium: e, inclusion: air
	if ( phi_a > 0.0 ) {
		HanaiBruggeman<T> hb;
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	}
	else
		return e;
}

//
// equation 5 de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_5(const T phi_c, // -> argile
											  const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f, // -> frequence
											  const Fluide<T>& fl) const
{
	std::complex<T> e_w, e_c;
	if ( f > 0.1 ) {
		e_w = ColeCole_w(f, fl.getSigma(), pPhys.temperature, fl.getSalinite())*iepsilon0;
		T tmp = 1./pPhys.fr_c;
		e_c = ColeCole(pPhys.Ks_c, pPhys.Ki_c, tmp, f,
					   pPhys.ccq_c, sigma_c*iepsilon0);
	}
	else {
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(fl.getSalinite());
		e_c = pPhys.Ks_c;
	}
	
	// medium: eau, inclusion: argile
	T phi = phi_w+phi_c;
	MaxwellGarnett<T> mg(35.0);
	std::complex<T> e = mg( e_c, e_w, phi_w/phi );
	// medium: sable, inclusion: e
	phi += phi_s;
	e = mg(e, std::complex<T>(pPhys.Km, 0.0), phi_s/phi);
	// medium: e, inclusion: air
	if ( phi_a > 0.0 ) {
		HanaiBruggeman<T> hb;
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	}
	else
		return e;
}

//
// equation 6 de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_6(const T phi_c, // -> argile
											  const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f // -> frequence
											  ) const
{
	std::complex<T> e_w, e_c;
	if ( f > 0.1 )
	{
		e_w = ColeCole_w(f, sw, pPhys.temperature,
                         Pride<T>::ppm2mol_l(pPhys.TDS)) *iepsilon0;
		T tmp = 1./pPhys.fr_c;
		e_c = ColeCole(pPhys.Ks_c, pPhys.Ki_c, tmp, f, pPhys.ccq_c,
					   sigma_c*iepsilon0);
	}
	else
	{
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(Pride<T>::ppm2mol_l(pPhys.TDS));
		e_c = pPhys.Ks_c;
	}
	
	// medium: argile, inclusion: eau
	T phi = phi_w+phi_c;
	HanaiBruggeman<T> hb;
	std::complex<T> e = hb( e_w, e_c, 1./3., phi_c/phi );
	// medium: e, inclusion: sable
	e = hb( std::complex<T>(pPhys.Km, 0.0), e, 1./3., phi/(phi+phi_s) );
	phi += phi_s;
	// medium: e, inclusion: air
	if ( phi_a > 0.0 )
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	else
		return e;
}

//
// equation 6 de Druchinin
template<typename T, template<typename>class GEN, typename COV>
std::complex<T> MilieuPoreux<T,GEN,COV>::eq_6(const T phi_c, // -> argile
											  const T phi_s, // -> sable
											  const T phi_w, // -> eau
											  const T phi_a, // -> air
											  const T f, // -> frequence
											  const Fluide<T>& fl) const
{
	std::complex<T> e_w, e_c;
	if ( f > 0.1 )
	{
		e_w = ColeCole_w(f, fl.getSigma(), pPhys.temperature, fl.getSalinite()) *iepsilon0;
		T tmp = 1./pPhys.fr_c;
		e_c = ColeCole(pPhys.Ks_c, pPhys.Ki_c, tmp, f, pPhys.ccq_c,
					   sigma_c*iepsilon0);
	}
	else
	{
		e_w = MalmbergMaryott( pPhys.temperature ) + corr_e_s(fl.getSalinite());
		e_c = pPhys.Ks_c;
	}
	
	// medium: argile, inclusion: eau
	T phi = phi_w+phi_c;
	HanaiBruggeman<T> hb;
	std::complex<T> e = hb( e_w, e_c, 1./3., phi_c/phi );
	// medium: e, inclusion: sable
	e = hb( std::complex<T>(pPhys.Km, 0.0), e, 1./3., phi/(phi+phi_s) );
	phi += phi_s;
	// medium: e, inclusion: air
	if ( phi_a > 0.0 )
		return hb(std::complex<T>(1.0, 0.0), e, 1./3., phi/(phi+phi_a) );
	else
		return e;
}

#endif
