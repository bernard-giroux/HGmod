//  Définition des constantes
//
//  Bernard Giroux
//  Ecole Polytechnique
//  29-11-2000

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

#ifndef __CONSTANTES_H__
#define __CONSTANTES_H__

#include <cmath>
#include <complex>

namespace Utils_bg {

  //  Pi
  const double pi = 4.0*atan(1.0);

  //  Permittivité du vide
  const double epsilon0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;        // farads/mètre
  const double eps0 = 8.851878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e+18;

  //  Perméabilité du vide
  const double mu0 = 4*pi*1e-7;             // Henry/mètre
  const double amu0 = 4*pi*1e-23;

  //  Vitesse de la lumière
  const double c0 = 299792458.0;            // mètre/seconde

  //  Charge electrique
  const double e = 1.6022e-19;              // Coulomb

  //  Constante de Boltzmann
  const double k = 1.38e-23;                // J/K

  //  Nombre d'Avogadro
  const double Na = 6.022e23;               // 1/mol
  
  const double small = 2.220446049250313e-16;
    
  const std::complex<double> i32 = std::pow(std::complex<double>(0.0, 1.0), 1.5);
};
#endif   // CONSTANTES_H
