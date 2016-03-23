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

#ifndef __UTILS_BG_H__
#define __UTILS_BG_H__

#include <complex>
#include <string>

#ifdef FFTW2
extern "C" {
#include "fftw.h"
}
#else
extern "C" {
#include "fftw3.h"
}
#endif

#include "enum_bg.h"


namespace Utils_bg {

  const int exit_fnf = 1;    // file not found
  const int exit_ni  = 2;    // non imlémenté
  const int exit_pff = 3;    // problème avec le format du fichier
  const int exit_pi  = 4;    // paramètres inconsistants
  const int exit_esf = 5;    // erreur de sortie fichier
  const int exit_mi  = 6;    // mémoire insuffisante

  struct uint_xyz {
	unsigned int x;
	unsigned int y;
	unsigned int z;
  };

  struct double_xyz {
	double       x;
	double       y;
	double       z;
  };

  struct uint_ijk {
	unsigned int i;
	unsigned int j;
	unsigned int k;
  };

  struct size_t_ijk {
	size_t       i;
	size_t       j;
	size_t       k;
  };

  struct double_ijk {
	double       i;
	double       j;
	double       k;
  };

  struct float_ijk_p {
	float*      i;
	float*      j;
	float*      k;
  };

  struct double_ijk_p {
	double*      i;
	double*      j;
	double*      k;
  };

  struct cplx_float_ijk_p {
	std::complex<float>*      i;
	std::complex<float>*      j;
	std::complex<float>*      k;
  };

  struct cplx_double_ijk_p {
	std::complex<double>*      i;
	std::complex<double>*      j;
	std::complex<double>*      k;
  };

  struct floatInt_xyz {
	float        x;
	float        y;
	float        z;
	int          i;
	int          j;
	int          k;
  };

  struct int_xyz {
	int          i;
	int          j;
	int          k;
  };

  struct float_xyz {
	float        x;
	float        y;
	float        z;
  };

  struct fftwf_complex_ijk_p {
	fftwf_complex* i;
	fftwf_complex* j;
	fftwf_complex* k;
  };

  struct fftw_complex_ijk_p {
	fftw_complex* i;
	fftw_complex* j;
	fftw_complex* k;
  };


  struct GMT_LUT {
	double z_low, z_high, i_dz;
	int rgb_low[3], rgb_high[3], rgb_diff[3];
	int anot;
	bool skip;
  };

};
#endif
