/*
 *  Fluide.h
 *  HGmod
 *
 *  Created by Bernard Giroux on 08-01-18.
 *  Copyright 2008 Bernard Giroux. All rights reserved.
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

#ifndef __FLUIDE_H__
#define __FLUIDE_H__

#include <cmath>
#include <complex>

#include "constantes.h"

template<typename T>
class Pride;

// fluide comportant 2 ions
template<typename X>
class Fluide {
public:
    // s: salinité [mol/l]
    // T: température [K]
    // z: valence
    // b: mobilité [(m/s)/N]
    // rho: densité [kg/m3]
    // eta: viscosité [cP]
    // epsilon: permittivité diélectrique [F/m]
    Fluide( X s_, X t, X rho_, X eta_, X epsilon_, X pH_, X z_[], X b_[])
    : s(s_), T(t), rho(rho_), eta(eta_), epsilon(epsilon_), pH(pH_), sigma(0),
    d(0)
    {
        z[0] = z_[0];
        z[1] = z_[1];
        ez[0] = Utils_bg::e*z[0];
        ez[1] = Utils_bg::e*z[1];
        b[0] = b_[0];
        b[1] = b_[1];
        beta[0] = fabs(ez[0])*b[0];
        beta[1] = fabs(ez[1])*b[1];
        B[0] = beta[0] + 2.0*epsilon*Utils_bg::k*T/(eta*ez[0]);
        B[1] = beta[1] + 2.0*epsilon*Utils_bg::k*T/(eta*ez[1]);
        D[0] = b[0]*Utils_bg::k*T;
        D[1] = b[1]*Utils_bg::k*T;
        update_tout();
    }
    
    Fluide(const Fluide<X>& fl)
    : s(fl.s), T(fl.T), rho(fl.rho), eta(fl.eta), epsilon(fl.epsilon), pH(fl.pH),
    sigma(0), d(0)
    {
        z[0] = fl.z[0];
        z[1] = fl.z[1];
        ez[0] = fl.ez[0];
        ez[1] = fl.ez[1];
        b[0] = fl.b[0];
        b[1] = fl.b[1];
        beta[0] = fl.beta[0];
        beta[1] = fl.beta[1];
        B[0] = fl.B[0];
        B[1] = fl.B[1];
        D[0] = fl.D[0];
        D[1] = fl.D[1];
        update_tout();
    }
    
    void setSalinite(const X s_) {
        s = s_;
        update_tout();
    }
    
    size_t getIndexContreIon() const { return ez[0] > ez[1] ? 0 : 1; }
    
    X getSalinite() const { return s; }
    
    X getSigma() const{ return sigma; }
    
    std::complex<X> getSigma(X w) const { return std::complex<X>(sigma, w*epsilon); }
    
    X getEpsilon() const { return epsilon; }
    
    std::complex<X> getEpsilon(X w) const { return std::complex<X>(epsilon, sigma/w); }

    const X* getCharge() const { return ez; }
    
    const X* getMobilite() const { return b; }
    
    const X* getMobiliteEffective() const { return B; }

    const X* getHittorf() const { return t_hf; }
    
    const X* getCoeffDiffusion() const { return D; }
    
    const X* getConcentration() const { return N; }
    
    X getEpaisseurDebye() const { return d; }
    
    X getTemperature() const { return T; }
    void setTemperature(const X t) {
        T = t;
        update_tout();
    }
    
    X getViscosite() const { return eta; }
    
    X getDensite() const { return rho; }
    
    // p. 4103 de Sen, 1987, J. Chem. Phys., vol. 87
    X getSurfaceCharge(X zeta) const {
        return 2.0 * d * N[getIndexContreIon()] * std::exp( 0.5*zeta );
    }
    
    X getSurfaceCharge() const {
        return getSurfaceCharge( Pride<X>::calculZeta( s ) );
    }

private:

    X s;          // salinité [mol/l]
    X T;          // temperature [K]
    X rho;        // masse volumique [kg/m3]
    X eta;        // viscosité [cP]
    X epsilon;    // permittivité diélectrique [F/m]
    X pH;
    X z[2];       // valence
	X ez[2];      // charge (C)
    X b[2];       // mobilité [(m/s)/N]
    X beta[2];    // mobilité [(m2/s)/V]
    X B[2];       // mobilité effective
    X N[2];       // concentration ionique (bulk) [1/m3]
    X t_hf[2];    // Hittorf Numbers
    X D[2];       // coefficient de diffusion
    X sigma;      // conductivité électrique [S/m]
    X d;          // epaisseur de Debye [m]
    
    // N: facteur de 1000 pour passer de 1/litre à 1/m3
    void update_N() { N[0] = N[1] = 1000.*s*Utils_bg::Na; }
    
    // Epaisseur de Debye, eq A-5 de Carcione
    void update_d()
    {
//        d = 1./std::sqrt((ez[0]*ez[0]*N[0]+ez[1]*ez[1]*N[1])/
//                         (epsilon*Utils_bg::k*T) );
        d = std::sqrt( (epsilon*Utils_bg::k*T)/
                      (ez[0]*ez[0]*N[0]+ez[1]*ez[1]*N[1])  );
    }
    
   void update_tout() {
        update_N();
        update_d();
        sigma = (ez[0]*ez[0]*b[0]*N[0]) + (ez[1]*ez[1]*b[1]*N[0]);
        t_hf[0] = ez[0] * beta[0] * N[0] / sigma;
        t_hf[1] = ez[1] * beta[1] * N[1] / sigma;
    }
    
};

#endif