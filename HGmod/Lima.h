/*
 *  Schurr.h
 *  HGmod
 *
 *  Created by Bernard Giroux on 08-01-19.
 *  Copyright 2008 École Polytechnique de Montréal. All rights reserved.
 *
 */
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
 *
 * @ARTICLE{lima92,
 *  author = {Olivar A. L. de Lima and Mukul M. Sharma},
 *  title = {A generalized {M}axwell-{W}agner theory for membrane polarization in shaly sands},
 *  journal = g,
 *  year = {1992},
 *  volume = {57},
 *  pages = {431--440},
 *  number = {3},
 *  doi = {10.1190/1.1443257}
 * }
 *
 */

#ifndef __LIMA_H__
#define __LIMA_H__

#include <complex>
#include <vector>

#include "constantes.h"
#include "Fluide.h"

template<typename T>
class Lima {
public:
	Lima(const Fluide<T>* fl) : fluide(fl) { }
	
	std::complex<T> sigma_schurr(const T a, const T omega);
	std::complex<T> sigma_fixman(const T a, const T omega);

    void setRho(T rho_b, T rho_d) { updateLambda(rho_b, rho_d); }
//	T s_model_sigma_ef(const T a, const T omega);
//	T s_model_epsilon_ef(const T a, const T omega);
//	
//	T d_model_sigma_ef(const T a, const T omega);
//	T d_model_epsilon_ef(const T a, const T omega);
	
private:
	const Fluide<T>* fluide;
    T lambda;
    T lambda0;
    
    void updateLambda(T rho_b, T rho_d)
    {
        T mu = Utils_bg::e * fluide->getMobilite()[ fluide->getIndexContreIon() ];
        lambda  = Utils_bg::e * rho_d * mu;
        lambda0 = Utils_bg::e * rho_b * mu;
    }
    
    T getTau(const T a)
    {
        return 0.5*a*a / fluide->getCoeffDiffusion()[ fluide->getIndexContreIon() ];
    }
};

template<typename T>
std::complex<T> Lima<T>::sigma_schurr(const T a, const T omega)
{
    T tau = getTau(a);
    T w2t2 = omega*omega*tau*tau;
    T re = 2*(lambda + (lambda0*w2t2)/(1.0+w2t2))/a;
    T im = -omega*(fluide->getEpsilon() + (2.0*lambda0*tau)/(a*(1.0+w2t2)));
    return std::complex<T>(re, im);
}

template<typename T>
std::complex<T> Lima<T>::sigma_fixman(const T a, const T omega)
{
    T tau = getTau(a);
    
    T cte1 = std::sqrt(omega*tau);
    std::complex<T> cte2 = 1.0 + std::complex<T>(1.0, -1.0) * cte1;
    std::complex<T> Y = -(cte2)/(2.0 * cte2 + std::complex<T>(0.0, 2.8*omega*tau));

    T beta1 = fluide->getSurfaceCharge();
    T delta1 = beta1 / (a * fluide->getConcentration()[ fluide->getIndexContreIon() ]);
    return delta1 * fluide->getSigma(omega) / (1.0 - delta1*Y);
}


#endif