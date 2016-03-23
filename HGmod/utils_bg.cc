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

#include <iostream>
#include <fstream>

#include "utils_bg_fct.h"
#include "constantes.h"

using namespace std;
using namespace Utils_bg;

// -----------------------------------------------------------------
// Utils_bg::drand(int &idum)
//
// retourne un nombre aléatoire entre 0.0 et 1.0 exclusivement,
// distribution uniforme
// -----------------------------------------------------------------
double Utils_bg::drand(int &idum) {
    const int IA=16807, IM=2147483647, IQ=127773, IR=2836, NTAB=32;
    const int NDIV=(1+(IM-1)/NTAB);
    const double EPS=3.0e-16, AM=1.0/IM,RNMX=(1.0-EPS);
    static int iy=0;
    static std::valarray<int> iv(NTAB);
    int j, k;
    double temp;
    
    if (idum <= 0 || !iy) {
        if (-idum < 1) idum=1;
        else idum = -idum;
        for (j=NTAB+7; j>=0; j--) {
            k = idum/IQ;
            idum = IA*(idum-k*IQ)-IR*k;
            if (idum < 0) idum += IM;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
    }
    k = idum/IQ;
    idum = IA*(idum-k*IQ)-IR*k;
    if (idum < 0) idum += IM;
    j = iy/NDIV;
    iy = iv[j];
    iv[j] = idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

// -----------------------------------------------------------------
// Utils_bg::randn(int &idum)
// -----------------------------------------------------------------
double Utils_bg::randn(int &idum) {
    static int iset=0;
    static double gset;
    double v1, v2, rsq, fac;
    if (idum < 0) iset = 0;
    if ( iset == 0 ) {
        do {
            v1 = 2.0*drand(idum)-1.0;
            v2 = 2.0*drand(idum)-1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    }
    else {
        iset = 0;
        return gset;
    }
}

double Utils_bg::bg_erf(double x)
{
    static int Nmax = 50;
    double tmp = 0.;
    for (int n=0; n<Nmax; ++n) {
        int s = (n%2) == 0 ? 1 : -1;
        tmp += s*pow(x, 2*n+1)/(bg_factoriel<double>(n)*(2*n+1));
    }
    return 2.*tmp/sqrt(pi);
}
