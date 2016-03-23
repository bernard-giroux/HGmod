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

#ifndef __UTILS_BG_FCT_H__
#define __UTILS_BG_FCT_H__

#include <valarray>
#ifdef MACSTL
#include <macstl/valarray.h>
#endif

#include "utils_bg.h"

namespace Utils_bg {

	double randn(int &idum);
	double drand(int &idum);

#ifdef MACSTL
	template<typename T>
	void fenetre_erf(stdext::valarray<T>&);
	template<typename T>
	void fenetre_erf(stdext::valarray<T>&, int);
	template<typename T>
	void fenetre_erf(stdext::valarray<T>&, double, double);
	template<typename T>
	void fenetre_erf(stdext::valarray<T>&, int, double, double);
#endif
	template<typename T>
	void fenetre_erf(std::valarray<T>&);
	template<typename T>
	void fenetre_erf(std::valarray<T>&, int);
	template<typename T>
	void fenetre_erf(std::valarray<T>&, double, double);
	template<typename T>
	void fenetre_erf(std::valarray<T>&, int, double, double);
	
	
	double bg_erf(double);
	
	template <typename T>
	T bg_factoriel(T);
	
	template<int N>
	class bg_factorial {
	public:
		enum { result = N * bg_factorial<N-1>::result };
	};
	template<>
	class bg_factorial<2> {
	public:
		enum { result = 2 };
	};
	
	template <typename T>
	T masseVolEau(T, T);
	template <typename T>
	T viscositeEau(const T);
	
	template<typename T>
	T Archie(const T sigma_w, const T phi, const T a, const T m);
	
};

template <typename T>
inline T Utils_bg::bg_factoriel(T n)
{
  return ( n <= 1 ) ? 1 : n * bg_factoriel(n-1);
}

template<typename T>
void Utils_bg::fenetre_erf(std::valarray<T>& d)
{
  // retourne erf de 0 a 1 (plutot que -1 a 1)
  double pas = 2*exp(1.0)/(d.size()-1);
  for (size_t n=0; n<d.size(); ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas - exp(1.0) );
  }
}

template<typename T>
void Utils_bg::fenetre_erf(std::valarray<T>& d, int N)
{
  // retourne erf de 0 a 1 (plutot que -1 a 1)
  double pas = 2*exp(1.0)/(N-1);
  for (int n=0; n<N; ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas - exp(1.0) );
  }
}

template<typename T>
void Utils_bg::fenetre_erf(std::valarray<T>& d, double min, double max)
{
  // retourne erf normalisé, de 0 a 1 (plutot que -1 a 1)
  double pas = (max-min)/(d.size()-1);
  for (size_t n=0; n<d.size(); ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas + min );
  }
}

template<typename T>
void Utils_bg::fenetre_erf(std::valarray<T>& d, int N,
						   double min, double max)
{
  // retourne erf normalisé, de 0 a 1 (plutot que -1 a 1)
  double pas = (max-min)/(N-1);
  for (int n=0; n<N; ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas + min );
  }
}

#ifdef MACSTL
template<typename T>
void Utils_bg::fenetre_erf(stdext::valarray<T>& d)
{
  // retourne erf de 0 a 1 (plutot que -1 a 1)
  double pas = 2*exp(1.0)/(d.size()-1);
  for (size_t n=0; n<d.size(); ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas - exp(1.0) );
  }
}

template<typename T>
void Utils_bg::fenetre_erf(stdext::valarray<T>& d, int N)
{
  // retourne erf de 0 a 1 (plutot que -1 a 1)
  double pas = 2*exp(1.0)/(N-1);
  for (size_t n=0; n<N; ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas - exp(1.0) );
  }
}

template<typename T>
void Utils_bg::fenetre_erf(stdext::valarray<T>& d, double min, double max)
{
  // retourne erf normalisé, de 0 a 1 (plutot que -1 a 1)
  double pas = (max-min)/(d.size()-1);
  for (size_t n=0; n<d.size(); ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas + min );
  }
}

template<typename T>
void Utils_bg::fenetre_erf(stdext::valarray<T>& d, int N,
						   double min, double max)
{
  // retourne erf normalisé, de 0 a 1 (plutot que -1 a 1)
  double pas = (max-min)/(N-1);
  for (size_t n=0; n<N; ++n) {
    d[n] = 0.5 + 0.5*bg_erf( n*pas + min );
  }
}
#endif

template<typename T>
T Utils_bg::masseVolEau(T TDS, T temp)
{
	// TDS en mg/l
	// teperature en C
	// retourne la masse volumique (density) en kg/m3
	T rho = 1000.0;
	rho *= 1.0-(temp+288.9414)/(508929.2*(temp+68.12963))*(temp-3.9863)*(temp-3.9863);
	TDS *= 0.001;
	T A = 0.824493 - 0.0040899*temp + 0.000076438*temp*temp-0.00000082467*temp*temp*temp + 
		0.0000000053675*temp*temp*temp*temp;
	T B =  -0.005724 + 0.00010227*temp - 0.0000016546*temp*temp;
	rho += A*TDS + B*std::pow(TDS,(3/2)) + 0.00048314*TDS*TDS;
	return rho;
}

template<typename X>
X Utils_bg::viscositeEau(const X T)
{
	// formules de Weast
	// T en Celsius
	// viscosité en cP
	return (T<20) ? 100.*exp(2.303*(1301./(998.333+8.1855*(T-20.)+0.00585*(T-20.)*(T-20.)) - 3.30233)) : 
	1.002*exp(2.303*(1.3272*(20.-T)-0.001053*(T-20.)*(T-20.))/(T+105.));
}

template<typename X>
X KozenyCarman(X phi, X Ss, X T, X facteur=2.0)
{
	// eq 2-39 de Schön
	
	// phi: porosité
	// Ss: surface spécifique [1/m]
	// T: tortuosité
	// facteur: cte liée à la géométrie
	// retourne la perméabilité en mD
	return 1.013249996628410e+15*phi*phi*phi / (facteur*(1.-phi)*(1.-phi)*Ss*Ss*T*T);
}

template<typename X>
X KozenyCarmanPerco(X phi, X phi_p, X Ss, X B=1.0)
{
	// eq 1 de Mavko97
	
	// phi: porosité
	// phi_p: porosité seuil de percolation
	// Ss: surface spécifique [1/m]
	// facteur: cte liée à la géométrie
	// retourne la perméabilité en mD
	phi -= phi_p;
	return 1.013249996628410e+15*B*phi*phi*phi / (Ss*Ss*(1.-phi)*(1.-phi));
}

template <typename T>
T Utils_bg::Archie(const T sigma_w, const T phi, const T a, const T m) {
	return 1./a * sigma_w * std::pow(phi, m);
}


#endif
