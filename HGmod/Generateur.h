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

#ifndef __GENERATEUR_H__
#define __GENERATEUR_H__

#include <complex>
#include <valarray>

#include <cstdio>
#include <unistd.h>         // pour getopt et gethostname

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "fftw3.h"
}

#include "utils_bg_fct.h"

#include "constantes.h"
#include "structHGmod.h"
#include "Messages.h"

extern bool verbose;
extern Messages msg;

template<typename A>
class Generateur
{
public:
	Generateur(const paramGrille& g, const A* ac, const int s=-1)
	: pg(g), covar(ac), idum1(s) {}
	
	virtual ~Generateur() {}
	
	virtual double z(double x1, double x2, double x3) const = 0;
	
	virtual double get_moyenne() const = 0;
	
	virtual void reinitialise() = 0;
		
	double get_variance() const { return covar->get_variance(); }
    
    double get_covariance(const std::valarray<double> x1, const std::valarray<double> x2) const
    {
        return covar->C(x1, x2);
    }
    
    const A* get_covariance_fct() const { return covar; }
	
	virtual void affiche_param(std::ostream& out) const
	{
		out << "   - " << msg.getString("Générateur: Classe de base") << "\n";
		covar->affiche_param( out );
	}
	
protected:
	paramGrille pg;
	const A *covar;
	mutable int idum1;
};

/*
 *
 * @article{leravalec00,
 *    author = {Micka\"ele Le Ravalec and Beno\^{\i}t Noetinger and Lin Y. Hu},
 *    journal = mg,
 *    number = {6},
 *    pages = {701--723},
 *    title = {The {FFT} Moving Average ({FFT-MA}) Generator: An Efficient Numerical Method for Generating and Conditioning Gaussian Simulations},
 *    volume = {32},
 *    year = {2000}
 * }
 * 
 */

template<typename A>
class FFTMA : public Generateur<A>
{
public:
	FFTMA(const paramGrille& g, const A* ac, const int s=-1);
	~FFTMA() {
		delete[] z_;
		delete[] S;
		delete[] Z;
		fftw_destroy_plan(pz);
		fftw_destroy_plan(pgz);
	}
	
	double z(double x1, double x2, double x3) const;
	
	void reinitialise();
	
	double get_moyenne() const { return gz.sum()/gz.size(); }
	
	void affiche_param(std::ostream& out) const
	{
		out << "   - " << msg.getString("Générateur: FFT-MA") << "\n";
		Generateur<A>::covar->affiche_param( out );
	}
	
private:
	size_t NNN, Ni, Nj, Nk;
	double* z_;
	std::complex<double>* S;
	std::complex<double>* Z;
	fftw_plan pz;
	fftw_plan pgz;
	
	std::valarray<double> gz;
	void initialise();
};

template<typename A>
FFTMA<A>::FFTMA(const paramGrille& g, const A* ac, const int s)
: Generateur<A>(g, ac, s)
{

	if ( verbose ) {
		std::cout << "\n" << msg.getString("Création du générateur FFT-MA") << "\n";
		std::cout.flush();
	}
	
    Utils_bg::double_ijk d = Generateur<A>::pg.d;
    
	Ni = Generateur<A>::pg.nn.i;
	Nj = Generateur<A>::pg.nn.j;
	Nk = Generateur<A>::pg.nn.k;
	
    int fac = 4;
	std::valarray<double> x(3);

	if (Ni>1) {
		x[0] = Ni*d.i;
		x[1] = 0.;
		x[2] = 0.;
		while ( Generateur<A>::covar->C(x) > 1.e-6 ) {
			Ni *= 2;
			x[0] = Ni*d.i;
		}
	}
	if (Nj>1) {
		x[0] = 0.;
		x[1] = Nj*d.j;
		x[2] = 0.;
		while ( Generateur<A>::covar->C(x) > 1.e-6 ) {
			Nj *= 2;
			x[1] = Nj*d.j;
		}
	}
    if (Nk>1) {
		x[0] = 0.;
		x[1] = 0.;
		x[2] = Nk*d.k;
		while ( Generateur<A>::covar->C(x) > 1.e-6 ) {
			Nk *= 2;
			x[2] = Nk*d.k;
		}
	}
    
    if ( Generateur<A>::covar->faireRotation() ) {
        fac *= 3;
        if (Ni>1) Ni = Ni > fac*size_t(Generateur<A>::covar->get_a(0)/d.i) ? Ni :
            fac*size_t(Generateur<A>::covar->get_a(0)/d.i);
        if (Nj>1) Nj = Nj > fac*size_t(Generateur<A>::covar->get_a(1)/d.j) ? Nj :
            fac*size_t(Generateur<A>::covar->get_a(1)/d.j);
        if (Nk>1) Nk = Nk > fac*size_t(Generateur<A>::covar->get_a(2)/d.k) ? Nk :
            fac*size_t(Generateur<A>::covar->get_a(2)/d.k);
        
        std::valarray<double> NN(3);
        NN[0] = Ni;    NN[1] = Nj;    NN[2] = Nk;
        Generateur<A>::covar->rotation(NN);
        if (Ni>1) {
            Ni = Ni > size_t(NN[0]) ? Ni : size_t(NN[0]);
        }
        if (Nj>1) {
            Nj = Nj > size_t(NN[1]) ? Nj : size_t(NN[1]);
        }
        if (Nk>1) {
            Nk = Nk > size_t(NN[2]) ? Nk : size_t(NN[2]);
        }
    }
    
    if (Ni>1) Ni *= 2;
    if (Nj>1) Nj *= 2;
    if (Nk>1) Nk *= 2;

    size_t Ni2 = Ni/2;
    size_t Nj2 = Nj/2;
    size_t Nk2 = Nk/2;

    if ( verbose ) {
        std::cout << "  " << msg.getString("Taille de la grille d'échantillonnage de la fct de covariance") << ":\n";
        std::cout << '\t' << Ni << 'x' << Nj << 'x' << Nk << '\n';
        std::cout.flush();
    }
    
	NNN = Ni*Nj*Nk;

	double *C = new double[NNN];
	z_ = new double[NNN];
    S = new std::complex<double>[Ni*Nj*(Nk2+1)];
	Z = new std::complex<double>[Ni*Nj*(Nk2+1)];
	
	
	size_t nc=0;
	for ( size_t ni=0; ni<Ni; ++ni ) {
		x[0] = d.i*(ni<=Ni2 ? ni : -(double)(Ni-ni));
		for ( size_t nj=0; nj<Nj; ++nj ) {
			x[1] = d.j*(nj<=Nj2 ? nj : -(double)(Nj-nj));
			for ( size_t nk=0; nk<Nk; ++nk ) {
				x[2] = d.k*(nk<=Nk2 ? nk : -(double)(Nk-nk));
				z_[nc] = Utils_bg::randn(Generateur<A>::idum1);
				C[nc++] = Generateur<A>::covar->C(x);
			}
		}
	}
	
//	std::ofstream fout;
//	fout.open("cov.dat");
//	for ( size_t ni=0, nc=0; ni<Ni; ++ni ) {
//		x[0] = d.i*(ni<=Ni2 ? ni : -(double)(Ni-ni));
//		for ( size_t nk=0; nk<Nk; ++nk, ++nc ) {
//			x[2] = d.k*(nk<=Nk2 ? nk : -(double)(Nk-nk));
//			fout << x[0] << " " << x[2] << " " << C[nc] << std::endl;
//		}
//	}
//	fout.close();
	
#ifdef _OPENMP
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initialisation de fftw avec %d threads ... ", omp_get_max_threads());
        fflush(stdout);
    }
    fftw_init_threads();
    fftw_plan_with_nthreads( omp_get_max_threads() );
#else
    if ( verbose >= 1 ) {
        fprintf(stdout, "Initialisation de fftw ... ");
        fflush(stdout);
    }	
#endif

	char hostname[100];
	std::string wisdomfile;
	std::string home = getenv( "HOME" );
	gethostname( hostname, 100 );
	wisdomfile = home+"/.fftw/"+hostname+".wisdom.HGmod";
	FILE *wfile;
	if ( (wfile = fopen(wisdomfile.c_str(), "r")) != NULL ) {
		fftw_import_wisdom_from_file( wfile );
		fclose( wfile );
	}
	fftw_plan pC = fftw_plan_dft_r2c_3d(Ni, Nj, Nk, C,
                                        reinterpret_cast<fftw_complex*>(S),
										FFTW_MEASURE);
	pz  = fftw_plan_dft_r2c_3d(Ni, Nj, Nk, z_,
                               reinterpret_cast<fftw_complex*>(Z),
							   FFTW_MEASURE);
	pgz = fftw_plan_dft_c2r_3d(Ni, Nj, Nk,
                               reinterpret_cast<fftw_complex*>(Z), z_,
							   FFTW_MEASURE);
	
	if ( (wfile = fopen(wisdomfile.c_str(), "w")) != NULL ) {
		fftw_export_wisdom_to_file( wfile );
		fclose( wfile );
	}
	
	fftw_execute(pC);
	fftw_execute(pz);
    
	for ( size_t n=0; n<Ni*Nj*(Nk2+1); ++n ) {
		Z[n] = sqrt( S[n] ) * Z[n];
	}
	
	fftw_execute(pgz);
	
	double iNNN = 1./(NNN);
	
    delete[] C;
	gz.resize( Generateur<A>::pg.nn.i*Generateur<A>::pg.nn.j*Generateur<A>::pg.nn.k );
    
	for ( size_t ni=0; ni<Generateur<A>::pg.nn.i; ++ni )
		for ( size_t nj=0; nj<Generateur<A>::pg.nn.j; ++nj )
			for ( size_t nk=0; nk<Generateur<A>::pg.nn.k; ++nk )
				gz[(ni*Generateur<A>::pg.nn.j+nj)*Generateur<A>::pg.nn.k+nk] =
					iNNN*(z_[(ni*Nj+nj)*Nk+nk]);
	
	fftw_destroy_plan(pC);
	
	if ( verbose ) std::cout << msg.getString("Création du générateur terminée") << std::endl;	
}

template<typename A>
double FFTMA<A>::z(double xi, double xj, double xk) const
{
	size_t ni = size_t( 0.0000001 + (xi-Generateur<A>::pg.min.i)/Generateur<A>::pg.d.i );
	size_t nj = size_t( 0.0000001 + (xj-Generateur<A>::pg.min.j)/Generateur<A>::pg.d.j );
	size_t nk = size_t( 0.0000001 + (xk-Generateur<A>::pg.min.k)/Generateur<A>::pg.d.k );
    double pepite = ( Generateur<A>::covar->contient_pepite() ) ?
        pepite = Generateur<A>::covar->get_et_pepite() * Utils_bg::randn(Generateur<A>::idum1) :
        0.0;
    return gz[(ni*Generateur<A>::pg.nn.j+nj)*Generateur<A>::pg.nn.k+nk] + pepite;
}

template<typename A>
void FFTMA<A>::reinitialise()
{
    if ( verbose ) {
		std::cout << "\n" << msg.getString("Réinitialisation du générateur FFT-MA") << " ... ";
		std::cout.flush();
	}
    
	for (size_t n=0; n<NNN; ++n) z_[n] = Utils_bg::randn(Generateur<A>::idum1);
	fftw_execute(pz);
	
	size_t Nk2 = Nk/2;
	for ( size_t n=0; n<Ni*Nj*(Nk2+1); ++n ) {
		Z[n] = sqrt( S[n] ) * Z[n];
	}
	
	fftw_execute(pgz);
	
	double iNNN = 1./(NNN);
	
	for ( size_t ni=0; ni<Generateur<A>::pg.nn.i; ++ni )
		for ( size_t nj=0; nj<Generateur<A>::pg.nn.j; ++nj )
			for ( size_t nk=0; nk<Generateur<A>::pg.nn.k; ++nk )
				gz[(ni*Generateur<A>::pg.nn.j+nj)*Generateur<A>::pg.nn.k+nk] =
					iNNN*(z_[(ni*Nj+nj)*Nk+nk]);
    
    if ( verbose ) std::cout << msg.getString("terminée") << "\n";
	
}

/*
 * Methode des bandes tournantes
 *
 *
 *  @article{tompson89,
 *     author = {Andrew F. B. Tompson and Rachid Ababou and Lynn W. Gelhar},
 *     journal = wrr,
 *     number = {10},
 *     pages = {2227--2243},
 *     title = {Implementation of the Three-Dimensional Turning Bands 
 *              Random Field Generator},
 *     volume = {25},
 *     year = {1989}
 *  }
 *
 */

template<typename A>
class TBM : public Generateur<A>
{
public:
	TBM(const paramGrille& g, const A* ac, const int s=-1, const int LL=100);
	~TBM() {
		delete[] fftout;
		delete[] dW;
		fftw_destroy_plan(p);
	}
	
	double z(double x1, double x2, double x3) const;
	void reinitialise();
	double get_moyenne() const;
	void affiche_param(std::ostream& out) const {
		out << "   - Générateur: Bandes tournantes\n";
		Generateur<A>::covar->affiche_param( out );
	}
	
private:
	size_t L,M;
	double iL;
	std::valarray<double> zeta;
	double idzeta;
	std::valarray<double> u;
	std::valarray<double> k;
	std::valarray<double> Z;
	double dk;
	std::complex<double> *dW;
	std::complex<double> *fftout;
	fftw_plan p;

	void calcul_zeta();  // a appeler avec calcul_?() (M y est défini)
	void calcul_u();
	void calcul_k();
	void calcul_Z();
	
};

template<typename A>
TBM<A>::TBM(const paramGrille& g, const A* ac, const int s, const int LL)
: Generateur<A>(g, ac,s), L(LL)
{
	if ( verbose ) {
		std::cout << "\n" << msg.getString("Création du générateur par bandes tournantes") << " ... ";
		std::cout.flush();
	}
	
	calcul_zeta();
	calcul_k();
	calcul_u();
	
	dW = new std::complex<double>[M];
	fftout = new std::complex<double>[M];
	Z.resize(L*M);
	
	char hostname[100];
	std::string wisdomfile;
	std::string home = getenv( "HOME" );
	gethostname( hostname, 100 );
	wisdomfile = home+"/.fftw/"+hostname+".wisdom.1d";
	FILE *wfile;
	if ( (wfile = fopen(wisdomfile.c_str(), "r")) != NULL ) {
		fftw_import_wisdom_from_file( wfile );
		fclose( wfile );
	}
	p = fftw_plan_dft_1d(M, reinterpret_cast<fftw_complex*>(dW),
						 reinterpret_cast<fftw_complex*>(fftout),
						 FFTW_BACKWARD, FFTW_MEASURE);
	if ( (wfile = fopen(wisdomfile.c_str(), "w")) != NULL ) {
		fftw_export_wisdom_to_file( wfile );
		fclose( wfile );
	}
	
	calcul_Z();
	iL = 1./sqrt(1.0*L);
	
	if ( verbose ) std::cout << msg.getString("terminée") << std::endl;
}

template<typename A>
void TBM<A>::calcul_zeta()
{
	
	// facteurs d'échelle pour satisfaire l'échantillonnage
	//   du processus aléatoire
	double dx_min_fac = 0.25;
	double x_max_fac = 5.0;
	double lambda_min_fac = 0.05;
	double lambda_max_fac = 50.0;
	
	double lambda_min = Generateur<A>::covar->get_a_min();
	double lambda_max = Generateur<A>::covar->get_a_max();
	double dzeta = Generateur<A>::pg.d.i < Generateur<A>::pg.d.j ? Generateur<A>::pg.d.i : Generateur<A>::pg.d.j;
	dzeta = dzeta < Generateur<A>::pg.d.k ? dzeta : Generateur<A>::pg.d.k;
	dzeta *= dx_min_fac;
	dzeta = dzeta < lambda_min_fac*lambda_min ? dzeta :
		lambda_min_fac*lambda_min;
	
	std::valarray<double> vtmp(3);
	vtmp[0] = fabs(Generateur<A>::pg.min.i) > fabs(Generateur<A>::pg.max.i) ? fabs(Generateur<A>::pg.min.i) : fabs(Generateur<A>::pg.max.i);
	vtmp[1] = fabs(Generateur<A>::pg.min.j) > fabs(Generateur<A>::pg.max.j) ? fabs(Generateur<A>::pg.min.j) : fabs(Generateur<A>::pg.max.j);
	vtmp[2] = fabs(Generateur<A>::pg.min.k) > fabs(Generateur<A>::pg.max.k) ? fabs(Generateur<A>::pg.min.k) : fabs(Generateur<A>::pg.max.k);
	vtmp *= vtmp;
	double zeta_max = x_max_fac*sqrt( vtmp.sum() );
	zeta_max = (zeta_max > lambda_max_fac * lambda_max ? zeta_max :
				lambda_max_fac * lambda_max);
	M = int(1.000001 + zeta_max/dzeta);
	zeta_max = (M-1)*dzeta;
	zeta.resize(2*M-1);
	M = zeta.size();
	for ( size_t n=0; n<zeta.size(); ++n ) zeta[n] = -zeta_max+n*dzeta;
	idzeta = 1./dzeta;
}

template<typename A>
void TBM<A>::calcul_k()
{
	double dzeta = zeta[1] - zeta[0];
	dk = 2*Utils_bg::pi/(M*dzeta);
	k.resize(M);
	for ( size_t n=0; n<M/2; ++n ) { k[n] = n; k[M-n-1] = n+1; }
	k[M/2] = M/2;
	k *= dk;
}


template<typename A>
void TBM<A>::calcul_u()
{
	if (u.size() != 3*L)
		u.resize(3*L);
	for ( size_t n=0; n<L; ++n ) {
		double phi = acos(1.-2.*Utils_bg::drand(Generateur<A>::idum1));        // eq. 29
		double theta = 2.*Utils_bg::pi*Utils_bg::drand(Generateur<A>::idum1);
		u[n*3+0] = sin(phi)*cos(theta);    // eq. 11
		u[n*3+1] = sin(phi)*sin(theta);
		u[n*3+2] = cos(phi);
	}
}


template<typename A>
void TBM<A>::calcul_Z()
{
	
	std::valarray<double> uu(3);
	//  double ipi = 1./Utils_bg::pi;
	for ( size_t l=0; l<L; ++l ) {
		uu[0] = u[l*3+0];  uu[1] = u[l*3+1];  uu[2] = u[l*3+2];
		for ( size_t m=0; m<M; ++m ) {
			// eq. 7
			double S1 = 2.*Utils_bg::pi*k[m]*k[m]*Generateur<A>::covar->S(k[m]*uu);
			// eq. 18
			double tmp = sqrt(2.0*S1*dk);
			double phi = 2.0*Utils_bg::pi*Utils_bg::drand(Generateur<A>::idum1);
			// eq. 16
			dW[m] = std::complex<double>(tmp*cos(phi), tmp*sin(phi));
		}
		fftw_execute(p);
		for ( size_t m=0; m<M; ++m ) {
			Z[l*M+m] = fftout[m].real();
		}
	}
}

template<typename A>
inline double TBM<A>::z(double x1, double x2, double x3) const
{
	double z_ = 0.;
	for (size_t l=0; l<L; ++l) {
		double zetaTmp = x1*u[l*3+0] + x2*u[l*3+1] + x3*u[l*3+2];
		z_ += Z[l*M+int( 0.00001+idzeta*(zetaTmp-zeta[0]) )];
	}
	z_ *= iL;
	return z_;
}

template<typename A>
void TBM<A>::reinitialise()
{
    if ( verbose ) {
		std::cout << "\n" << msg.getString("Réinitialisation du générateur par bandes tournantes") << " ... ";
		std::cout.flush();
	}
    
	calcul_u();
	calcul_Z();
    
    if ( verbose ) std::cout << msg.getString("terminée") << std::endl;
}

template<typename A>
double TBM<A>::get_moyenne() const
{
	double moy = 0.;
	for ( size_t ni=0; ni<Generateur<A>::pg.nn.i; ++ni ) {
		double x1 = Generateur<A>::pg.min.i + ni*Generateur<A>::pg.d.i;
		for ( size_t nj=0; nj<Generateur<A>::pg.nn.j; ++nj ) {
			double x2 = Generateur<A>::pg.min.j + nj*Generateur<A>::pg.d.j;
			for ( size_t nk=0; nk<Generateur<A>::pg.nn.k; ++nk ) {
				double x3 = Generateur<A>::pg.min.k + nk*Generateur<A>::pg.d.k;
				moy += this->z(x1, x2, x3);
			}
		}
	}
	return moy/(Generateur<A>::pg.nn.i*Generateur<A>::pg.nn.j*Generateur<A>::pg.nn.k);
}



#endif
