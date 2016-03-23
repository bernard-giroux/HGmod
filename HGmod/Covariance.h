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

#ifndef __COVARIANCE_H__
#define __COVARIANCE_H__

#include <iostream>
#include <valarray>

#include "constantes.h"

#include "structHGmod.h"
#include "Messages.h"

extern bool verbose;
extern Messages msg;

class Covariance
{
public:
    Covariance(const paramCovariance& p) : a(p.a), ia(1./(p.a+Utils_bg::small)),
    theta(p.theta), palier(p.palier),
    cangz(0.0), sangz(0.0), cangy(0.0),
    sangy(0.0), cangx(0.0), sangx(0.0), a2(0.0),
    faireRotation_(false), type(p.typeCov)
	{
        std::valarray<double> aa = a*a;
        a2 = aa.sum();
		if ( !(theta[0]==0.0 && theta[1]==0.0 && theta[2]==0.0) ) {
			faireRotation_ = true;
			cangz = std::cos(p.theta[2]/180.0*Utils_bg::pi);
			sangz =	std::sin(p.theta[2]/180.0*Utils_bg::pi);
			cangy = std::cos(p.theta[1]/180.0*Utils_bg::pi);
			sangy = std::sin(p.theta[1]/180.0*Utils_bg::pi);
			cangx = std::cos(p.theta[0]/180.0*Utils_bg::pi);
			sangx = std::sin(p.theta[0]/180.0*Utils_bg::pi);
		}
	}

	virtual ~Covariance() {}
	
	double get_a(int n) const { return a[n]; }
	double get_a_min() const { return a.min(); }
	double get_a_max() const { return a.max(); }
	double get_palier() const { return palier; }
    Utils_bg::typeCovariance get_type() const { return type; }
    
	bool faireRotation() const { return faireRotation_; }

	virtual double C(std::valarray<double> x) const = 0;
	virtual double C(std::valarray<double> x1, std::valarray<double> x2) const = 0;
	virtual double S(std::valarray<double> k) const = 0;
	virtual void affiche_param(std::ostream& out) const {
		out << "                  " << msg.getString("palier") << ": " << palier;
		out << "\n                  " << msg.getString("portées") << ":";
        for ( unsigned int n=0; n<a.size(); ++n ) out << " " << a[n];
		out << "\n                  " << msg.getString("θ anisotropie") << ":";
        for ( unsigned int n=0; n<a.size(); ++n ) out << " " << theta[n];
        out << '\n';
    }

	void rotation(std::valarray<double>& v) const
	{
		if ( faireRotation_ ) {
			double tmp0 = v[0];
			double tmp1 = v[1];
			double tmp2 = v[2];
			v[0] = tmp0*cangz*cangy + tmp1*sangz*cangy - tmp2*sangy;
		
			v[1] = tmp0*(-sangz*cangx+cangz*sangy*sangx);
			v[1] += tmp1*(cangz*cangx+sangz*sangy*sangx);
			v[0] += tmp2*cangy*sangx;
		
			v[2] = tmp0*(sangz*sangx+cangz*sangy*cangx);
			v[2] += tmp1*(-cangz*sangx+sangz*sangy*cangx);
			v[2] += tmp2*cangy*cangx;
		}
	}	
	
protected:
	std::valarray<double> a;     // portées
	std::valarray<double> ia;
	std::valarray<double> theta;
	double palier;
	double cangz;
	double sangz;
	double cangy;
	double sangy;
	double cangx;
	double sangx;
	double a2;
	bool faireRotation_;
    Utils_bg::typeCovariance type;

	// calcul le lag normalisé par la portée
	double calcul_h(std::valarray<double>& x) const {
		rotation(x);
        x *= ia;
		x *= x;    
		return sqrt( x.sum() );
	}
	double calcul_h(std::valarray<double>& x1, std::valarray<double>& x2) const {
		x1 -= x2;		rotation(x1);
		x1 *= ia;       x1 *= x1;
		return sqrt( x1.sum() );
	}
	double calcul_h2(std::valarray<double>& x) const {
		rotation(x);
        x *= ia;
		x *= x;    
		return x.sum();
	}
	double calcul_h2(std::valarray<double>& x1, std::valarray<double>& x2) const {
		x1 -= x2;		rotation(x1);
		x1 *= ia;       x1 *= x1;
		return x1.sum();
	}
	
	template<typename T>
		T prod(const std::valarray<T>& v) const
	{
		T t=v[0]; for (size_t n=1; n<v.size(); ++n) t*=v[n]; return t;
	}

	template<typename T>
		T norm2(std::valarray<T> v) const
	{
		v *= v;
		return sqrt( v.sum() );
	}
};

class Regional
{
public:
    Regional() : cov(), pepite(-1.0), et_pepite(0.0) {}
    
    ~Regional() {
        for (size_t n=0; n<cov.size(); ++n) delete cov[n];
    }

	void remplire(const std::vector<paramCovariance>&);

	void push_back(Covariance* c) {
        if ( c->get_type() == Utils_bg::PEPITE ) {
            pepite = c->get_palier();
			et_pepite = sqrt(pepite);
		}
		else
			cov.push_back(c);
    }
	
	double get_a(int n) const { 
        double val = cov[0]->get_a(n);
        for (size_t nn=1; nn<cov.size(); ++nn)
            val = val > cov[nn]->get_a(n) ? val : cov[nn]->get_a(n);;
        return val;
    }
    
    double get_a_min() const { 
        double val = cov[0]->get_a_min();
        for (size_t n=1; n<cov.size(); ++n)
            val = val < cov[n]->get_a_min() ? val : cov[n]->get_a_min();
        return val;
    }
    
	double get_a_max() const {
        double val = cov[0]->get_a_max();
        for (size_t n=1; n<cov.size(); ++n)
            val = val > cov[n]->get_a_max() ? val : cov[n]->get_a_max();
        return val;
    }
    
    bool contient_pepite() const { return pepite > 0.0; }
    
    double get_palier_pepite() const { return pepite; }
    double get_et_pepite() const { return et_pepite; }

    bool faireRotation() const {
        bool val = cov[0]->faireRotation();
        for (size_t n=1; n<cov.size(); ++n)
            val = val && cov[n]->faireRotation();
        return val;
    }

	double C(std::valarray<double> x) const {
		double val = cov[0]->C(x);
		for (size_t n=1; n<cov.size(); ++n)
                val += cov[n]->C(x);
		return val;
	}
	
	double C(std::valarray<double> x1, std::valarray<double> x2) const {
		double val = cov[0]->C(x1, x2);
		for (size_t n=1; n<cov.size(); ++n)
            val += cov[n]->C(x1, x2);
		return val;
	}
	
	double S(std::valarray<double> k) const {
		double val = cov[0]->S(k);
		for (size_t n=1; n<cov.size(); ++n)
                val += cov[n]->S(k);
		return val;
	}
		
	void affiche_param(std::ostream& out) const {
		out << msg.getString("Modèle régional de covariance comportant") << ":\n";
		for (size_t n=0; n<cov.size(); ++n) cov[n]->affiche_param(out);
	}
    
    // on retourne la longueur max pour chaque dimension, après rotation
    template<typename T>
    std::valarray<T> rotation(const std::valarray<T>& v) const {
        std::valarray<T> tmp = v;
		cov[0]->rotation(tmp);
        for (size_t n=1; n<cov.size(); ++n) {
            std::valarray<T> tmp2 = v;
			cov[n]->rotation(tmp2);
            for (size_t nn=0; nn<tmp2.size(); ++nn)
                tmp[nn] = tmp[nn] > tmp2[nn] ? tmp[nn] : tmp2[nn];
        }
        return tmp;
    }
    
private:
	std::vector<Covariance*> cov;
    double pepite, et_pepite;
};

class Pepite : public Covariance
{
public:
    Pepite(const paramCovariance& p) : Covariance(p) {
        type = Utils_bg::PEPITE;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance effet pépite") << '\n';
	}
    
    double C(std::valarray<double> x) const
    {
        return norm2(x) < Utils_bg::small ? palier : 0.0;
	}
    
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
        x1 -= x2;
        return norm2(x1) < Utils_bg::small ? palier : 0.0;
	}
    
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
        std::cerr << "Erreur: Densité spectrale du modèle de covariance gaussien non implémentée.\n";
        std::cerr << "Utilisez le générateur FFT-MA pour utiliser un modèle de covariance gaussien.\n";
        abort();
        return 0.0;
    }
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Effet pépite") << '\n';
		out << "                  " << msg.getString("palier") << ": " << palier << '\n';
    }
};

class VonKarman : public Covariance
{
public:
    VonKarman(const paramCovariance& p) : Covariance(p), v(p.nbHurst),
    E(a.size()), Snum(0.0), Ccte(0.0)
    {
        type = Utils_bg::VON_KARMAN;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance de von Karman") << '\n' ;
        double gamma_vE = gamma( v+0.5*E );
        double gamma_v = gamma( v );
        double cte1 = palier * std::pow(2.*sqrt(Utils_bg::pi), E)*prod(a);
        Snum = cte1*gamma_vE/gamma_v /(8.*Utils_bg::pi*Utils_bg::pi*Utils_bg::pi);
        Ccte = palier / (pow(2.0,v-1)*gamma_v);
    }
    
    // Covariance
    double C(std::valarray<double> x) const
    {
		double h = calcul_h(x);
		return h==0.0 ? palier : Ccte * pow( h, v ) * Kv( h );
    }
	
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
		double h = calcul_h(x1, x2);
        return h==0.0 ? palier : Ccte * pow( h, v ) * Kv( h );
    }
    
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
		rotation(k);
        k *= a;    k *= k;
        return Snum/pow(1. + k.sum(), v+0.5*E);
    }
    
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Modèle Von Karman") << '\n';
        Covariance::affiche_param(out);
        out << "\n                  " << msg.getString("Nb de Hurst") << ": " << v << '\n';
    }
    
private:
    
    double v;    // Hurst number
    double E;    // dimension Euclidienne
    
    double Snum, Ccte;
    
    double gammln(const double xx)
    {
        // Returns the value ln[ (xx)] for xx > 0.
        // Internal arithmetic will be done in double precision, a nicety that
        // you can omit if  five-figure accuracy is good enough.
        static const double cof[6] = {76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5395239384953e-5};
        
        double y=xx;
        double tmp=xx+5.5;
        tmp -= (xx+0.5)*std::log(tmp);
        double ser=1.000000000190015;
        for (int j=0; j<6; j++) ser += cof[j]/++y;
        return -tmp+std::log(2.5066282746310005*ser/xx);
    }
    double gamma(const double xx) { return std::exp( gammln(xx) ); }
    
    double Kv( const double ) const;
    void beschb(const double x, double &gam1, double &gam2,
                double &gampl, double &gammi) const;
    double chebev(const double a, const double b, std::valarray<double>& c,
                  const double x) const;
    
};


class Exponentiel : public Covariance
{
public:
    Exponentiel(const paramCovariance& p) : Covariance(p), Snum(0.0)
    {
        type = Utils_bg::EXPONENTIEL;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance exponentiel") << '\n';
        Snum = palier*prod(a)/(Utils_bg::pi*Utils_bg::pi);
    }
    
    double C(std::valarray<double> x) const
    {
		return palier*exp( -calcul_h(x) );
    }
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
		return palier*exp( -calcul_h(x1, x2) );
    }
	
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
		rotation(k);
        k *= a;    k *= k;    double tmp = 1. + k.sum();
        return Snum/(tmp*tmp);
    }
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Modèle Exponentiel") << '\n';
        Covariance::affiche_param(out);
    }
private:
    double Snum;
};


class Spherique : public Covariance
{
public:
    Spherique(const paramCovariance& p) : Covariance(p) 
    {
        type = Utils_bg::SPHERIQUE;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance sphérique") << '\n';
	}
    
    double C(std::valarray<double> x) const
    {
        double h = calcul_h(x);
        return palier*( h > 1.0 ? 0.0 : 1.0 - 1.5*h + 0.5*h*h*h);
	}
    
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
        double h = calcul_h(x1, x2);
        return palier*( h > 1.0 ? 0.0 : 1.0 - 1.5*h + 0.5*h*h*h);
	}
    
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
        std::cerr << "Erreur: Densité spectrale de la fct de covariance sphérique non implémentée.\n";
        std::cerr << "Utilisez le générateur FFT-MA pour utiliser un modèle de covariance sphérique.\n";
        abort();
        return 0.0;
    }
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Modèle Sphérique") << '\n';
        Covariance::affiche_param(out);
    }
};


class Gaussien : public Covariance
{
public:
    Gaussien(const paramCovariance& p) : Covariance(p), Snum(0.0)
    {
        type = Utils_bg::GAUSSIEN;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance gaussienne") << '\n';
		double tmp = sqrt(Utils_bg::pi);
		Snum = palier*prod(a) / (8. * tmp*tmp*tmp);
	}
    double C(std::valarray<double> x) const
    {
        return palier*exp( -calcul_h2(x) );
	}
    
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
        return palier*exp( -calcul_h2(x1, x2) );
	}
    
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
        // À vérifier...
		rotation(k);
        k *= a;    k *= k;
		return Snum * exp( -0.25 * k.sum() );
    }
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Modèle Gaussien") << '\n';
        Covariance::affiche_param(out);
    }
private:
	double Snum;
};

class Hyperbolique : public Covariance
{
public:
    Hyperbolique(const paramCovariance& p) : Covariance(p) 
    {
        type = Utils_bg::HYPERBOLIQUE;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance hyperbolique") << '\n';
	}
    
    double C(std::valarray<double> x) const
    {
        return palier/(1.+calcul_h(x));
	}
    
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
        return palier/(1.+calcul_h(x1, x2));
	}
    
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
        std::cerr << "Erreur: Densité spectrale de la fct de covariance hyperbolique non implémentée.\n";
        std::cerr << "Utilisez le générateur FFT-MA pour utiliser un modèle de covariance hyperbolique.\n";
        abort();
        return 0.0;
    }
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Modèle Hyperbolique") << '\n';
        Covariance::affiche_param(out);
    }
};

class Stable : public Covariance
{
public:
    Stable(const paramCovariance& p) : Covariance(p) 
    {
        type = Utils_bg::STABLE;
        if ( verbose )
            std::cout << '\n' << msg.getString("Création de la fonction de covariance Stable (α = 1/2)") << '\n';
	}
    
    double C(std::valarray<double> x) const
    {
        return palier*exp( -sqrt(calcul_h(x)) );
	}
    
    double C(std::valarray<double> x1, std::valarray<double> x2) const
    {
        return palier*exp( -sqrt(calcul_h(x1, x2)) );
	}
    
    // Densité Spectrale
    double S(std::valarray<double> k) const
    {
        std::cerr << "Erreur: Densité spectrale de la fct de covariance Stable non implémentée.\n";
        std::cerr << "Utilisez le générateur FFT-MA pour utiliser un modèle de covariance Stable.\n";
        abort();
        return 0.0;
    }
    void affiche_param(std::ostream& out) const
    {
        out << "   - " << msg.getString("Covariance") << " - " << msg.getString("Modèle Stable (α = 1/2)") << '\n';
        Covariance::affiche_param(out);
    }
};


#endif
 