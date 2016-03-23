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
#include <limits>

#include "Covariance.h"

extern Messages msg;

using namespace std;

void Regional::remplire(const std::vector<paramCovariance>& pCov) {
	for (size_t n=0; n<pCov.size(); ++n) {
		switch ( pCov[n].typeCov ) {
			case Utils_bg::EXPONENTIEL :
				cov.push_back( new Exponentiel( pCov[n] ) );
				break;
				
			case Utils_bg::VON_KARMAN :
				cov.push_back( new VonKarman( pCov[n] ) );
				break;
				
			case Utils_bg::SPHERIQUE :
				cov.push_back( new Spherique( pCov[n] ) );
				break;
				
			case Utils_bg::GAUSSIEN :
				cov.push_back( new Gaussien( pCov[n] ) );
				break;
				
			case Utils_bg::PEPITE :
                pepite = pCov[n].palier;
				et_pepite = sqrt(pepite);
				break;
				
			case Utils_bg::HYPERBOLIQUE :
				cov.push_back( new Hyperbolique( pCov[n] ) );
				break;
                
			case Utils_bg::STABLE :
				cov.push_back( new Stable( pCov[n] ) );
				break;
                
			default:
				cerr << msg.getString("Erreur, modèle de covariance non implémenté") << endl;
				abort();
		}
	}
}


double VonKarman::Kv(const double x) const {
    const int MAXIT = 10000;
    const double EPS = numeric_limits<double>::epsilon();
    const double FPMIN = numeric_limits<double>::min()/EPS;
    const double XMIN = 2.0, PI = 4.0*atan(1.0);
    
    if ( x <= 0.0 || v < 0.0 ) {
        cerr << "Erreur Kv: x = " << x << ", v = " << v << endl;
        exit(1);
    }
    int nl = int(v+0.5);
    double xmu = v-nl;
    double xmu2 = xmu*xmu;
    double xi = 1.0/x;
    double xi2 = 2.0*xi;
    double h = v*xi;
    if ( h < FPMIN ) h = FPMIN;
    double b = xi2*v;
    double d = 0.0;
    double c = h;
    int i;
    for ( i=0; i<MAXIT; ++i) {
        b += xi2;
        d = 1.0/(b+d);
        c = b+1.0/c;
        double del = c*d;
        h *= del;
        if ( fabs(del-1.0) <= EPS ) break;
    }
    if ( i >= MAXIT ) {
        cerr << "Erreur Kv: x trop grand" << endl;
        exit(1);
    }
    double fact = v*xi;
    double rkmu;
    double rk1;
    if ( x < XMIN ) {
        double x2 = 0.5*x;
        double pimu = PI*xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
        d = -log(x2);
        double e = xmu*d;
        double fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
        double gam1, gam2, gampl, gammi;
        beschb(xmu, gam1, gam2, gampl, gammi);
        double ff = fact*(gam1*cosh(e) + gam2*fact2*d);
        double sum = ff;
        e = exp(e);
        double p = 0.5*e/gampl;
        double q = 0.5/(e*gammi);
        c = 1.0;
        d = x2*x2;
        double sum1 = p;
        for ( i=1; i<=MAXIT; ++i ) {
            ff = (i*ff+p+q)/(i*i-xmu2);
            c *= (d/i);
            p /= (i-xmu);
            q /= (i+xmu);
            double del = c*ff;
            sum += del;
            double deli = c*(p-i*ff);
            sum1 += deli;
            if ( fabs(del) < fabs(sum)*EPS ) break;
        }
        if ( i > MAXIT ) {
            cerr << "Erreur Kv: convergence non atteinte" << endl;
            exit(1);
        }
        rkmu = sum;
        rk1 = sum1*xi2;
    } else {
        b = 2.0*(1.0+x);
        d = 1.0/b;
        double delh = h = d;
        double q1 = 0.0;
        double q2 = 1.0;
        double a1 = 0.25-xmu2;
        double q = c = a1;
        double a = -a1;
        double s = 1.0+q*delh;
        for ( i=1; i<MAXIT; ++i ) {
            a -= 2*i;
            c = -a*c/(i+1.0);
            double qnew = (q1-b*q2)/a;
            q1 = q2;
            q2 = qnew;
            q += c*qnew;
            b += 2.0;
            d = 1.0/(b+a*d);
            delh = (b*d-1.0)*delh;
            h += delh;
            double dels = q*delh;
            s += dels;
            if ( fabs(dels/s) <= EPS ) break;
        }
        if ( i >= MAXIT ) {
            cerr << "Erreur Kv: convergence non atteinte dans cf2" << endl;
            exit(1);
        }
        h *= a1;
        rkmu = sqrt(PI/(2.0*x))*exp(-x)/s;
        rk1 = rkmu*(xmu+x+0.5-h)*xi;
    }
    for ( i=1; i<=nl; ++i ) {
        double rktemp = (xmu+i)*xi2*rk1+rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    double rk = rkmu;
    //   rkp = v*xi*rkmu-rk1;
    return rk;
}

void VonKarman::beschb(const double x, double &gam1, double &gam2,
					   double &gampl, double &gammi) const
{
    static const double c1_d[7] = {
        -1.142022680371168e0,6.5165112670737e-3,
        3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
        3.67795e-11,-1.356e-13};
    static const double c2_d[8] = {
        1.843740587300905e0,-7.68528408447867e-2,
        1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
        2.423096e-10,-1.702e-13,-1.49e-15};
    static valarray<double> c1(c1_d,7), c2(c2_d,8);
    
    double xx = 8.0*x*x-1.0;
    gam1 = chebev(-1.0, 1.0, c1, xx);
    gam2 = chebev(-1.0, 1.0, c2, xx);
    gampl = gam2-x*gam1;
    gammi = gam2+x*gam1;
}

double VonKarman::chebev(const double a, const double b,
						 valarray<double>& c, const double x) const
{
    if ( (x-a)*(x-b) > 0.0 ) {
        cerr << "Erreur chebev: x en dehors de l'intervalle" << endl;
        exit(1);
    }
    double d=0.0, dd=0.0;
    double y = (2.0*x-a-b)/(b-a);
    double y2 = 2.0*y;
    double sv;
    for ( size_t j=c.size()-1; j>0; --j ) {
        sv = d;
        d = y2*d-dd+c[j];
        dd=sv;
    }
    return y*d-dd+0.5*c[0];
}

