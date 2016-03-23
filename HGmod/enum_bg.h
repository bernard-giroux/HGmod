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

#ifndef __ENUM_BG_H__
#define __ENUM_BG_H__

namespace Utils_bg {

  enum composante            { I_=0, J_=1, K_=2 };
  enum composantesChampTag   { Ex=0, Ey=1, Ez=2, Hx=3, Hy=4, Hz=5 };
  enum composantes           { IJK, IK, J };
  enum dimensionnalite       { _1D, _2D, _3D };
  enum fenetre               { COSH, BLACKMAN, HANNING, CK };
  enum fftLib                { FFTW, TEMPERTON };
  enum modeTransverse        { TE, TM, TEM }; // TEM -> 3D
  enum modeleConductivite    { PRIDE, REVIL };
  enum modeleRelaxation      { DEBYE };
  enum operateurDifferentiel { PSEUDOSPECTRAL, DIFFERENCES_FINIES };
  enum plan                  { XY, XZ, YZ };
  enum precision             { FLOAT, DOUBLE };
  enum typeChamp             { ELECTRIQUE, MAGNETIQUE, TMP };
  enum typeCoordonnees       { CARTESIENNES, CYLINDRIQUES };
  enum typeCovariance        { EXPONENTIEL, VON_KARMAN, SPHERIQUE, GAUSSIEN, PEPITE, HYPERBOLIQUE, STABLE };
  enum typeGenerateur        { TBM, FFTMA, DS };
  enum typeMaille            { REGULIERE, QUINCONCE };
  enum typePas               { UNIFORME, NON_UNIFORME };
  enum typeSource            { PONCTUELLE, LINEAIRE, PLANE };
  enum typeWavelet           { ASYM, SYM, RICKER, GAUSS, CHEN };
};

#endif
