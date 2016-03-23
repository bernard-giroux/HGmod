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

#include "Grille.h"

using namespace std;
using namespace Utils_bg;

//bool Grille::verbose=false;
//bool Grille::phi=false;
//bool Grille::Cl=false;

extern bool verbose;
extern Messages msg;

Grille::Grille(const paramGrille& pg, paramSortie& ps2) : p(pg), ps(ps2)
{
	
	if ( verbose ) cout << '\n' << msg.getString("Préparation de la grille") << endl;
	
	minm = p.min;
	maxm = p.max;
	p.nn.i = static_cast<size_t>( 1.00000001 + (maxm.i-minm.i) / p.d.i );
	p.nn.j = static_cast<size_t>( 1.00000001 + (maxm.j-minm.j) / p.d.j );
	p.nn.k = static_cast<size_t>( 1.00000001 + (maxm.k-minm.k) / p.d.k );
	
	if ( verbose ) {
		cout << " - " << msg.getString("La grille est de ") << p.nn.i << " x " << p.nn.j 
		     << " x " << p.nn.k << " " << msg.getString("noeuds") << ".\n";
		cout << " - " << msg.getString("Pas: ") << p.d.i<< " x " << p.d.j << " x " << p.d.k << '\n';
		cout << " - " << msg.getString("Domaine modélisé") << "   X: " << p.min.i << " - " << p.max.i
             << "\n                      Y: " << p.min.j << " - " << p.max.j
             << "\n                      Z: " << p.min.k << " - " << p.max.k
             << '\n';
	}
	
	size_t nn = p.nn.i*p.nn.j*p.nn.k;
    
    if (ps.sauvePorosite)       porosite.resize(nn);
	if (ps.sauvePermittivite)   permit.resize(nn);
    if (ps.sauveConductivite)   sigma.resize(nn);
    if (ps.sauvePermeaMag)      permea.resize(nn);
    if (ps.sauvePermeabilite)   permeabilite.resize(nn);
    if (ps.sauveFractionArgile) cl.resize(nn);
    if (ps.sauveChargeabilite)  chargeabilite.resize(nn);
    if (ps.sauveConductiviteComplexe) {
        sigmaCplx.resize(nn);
        for (size_t n=0; n<sigmaCplx.size(); ++n) sigmaCplx[n].resize( ps.nf );
    }
    if (ps.sauveVp) Vp.resize(nn);
	saliniteVariable = false;
    simulerPorosite = true;
}

void Grille::sauveGrille(const std::string fichier)
{
	if ( verbose ) {
		cout << '\n' << msg.getString("Sauvegarde de la grille") << " ... ";
		cout.flush();
	}
    ofstream fout((fichier+"_"+msg.getString("grille")+".dat").c_str());
    fout << "% x\ty\tz\n";
    for ( size_t ni=0, I=0; ni<p.nn.i; ++ni )
    {
        double x = minm.i + ni*p.d.i;
        for ( size_t nj=0; nj<p.nn.j; ++nj )
        {
            double y = minm.j + nj*p.d.j;
            for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
            {
                double z = minm.k + nk*p.d.k;
                fout << x << '\t' << y << '\t' << z << '\n';
            }
        }
    }
    fout.close();
	if ( verbose ) cout << msg.getString("terminée") << '\n';
}

void Grille::sauveSalinite(const std::string fichier)
{
	if ( verbose ) {
		cout << '\n' << msg.getString("Sauvegarde du modèle de salinité") << " ... ";
		cout.flush();
	}
    ofstream fout((fichier+"_"+msg.getString("salinité")+".dat").c_str());
    fout << msg.getString("% salinité") << '\n';
    for ( size_t n=0; n<salinite.size(); ++n )
        fout << salinite[n] << '\n';
    fout.close();
	if ( verbose ) cout << msg.getString("terminée") << '\n';
}


void Grille::conditionnePorosite(const paramConditionnement& pc)
{
    if ( verbose )
        std::cout << '\n' << msg.getString("Conditionnement des porosités simulées (krigeage simple)");
    
    size_t N = pc.coord.size();
	
    boost::numeric::ublas::vector<double> phi(N);
    for (size_t n=0; n<N; n++)
        phi[n] = pc.data[n];
    
    for ( size_t n=0; n<N; n++ )
    {
        size_t ni = static_cast<size_t>((0.0000001+pc.coord[n].i-p.min.i)/p.d.i);
        size_t nj = static_cast<size_t>((0.0000001+pc.coord[n].j-p.min.j)/p.d.j);
        size_t nk = static_cast<size_t>((0.0000001+pc.coord[n].k-p.min.k)/p.d.k);
        size_t ind = (p.nn.j*nk + nj)*p.nn.i + ni;
        
        // différence entre les données et valeurs simulées
        phi[n] -= porosite[ind];
    }
    
    if ( verbose ) std::cout << "\n   " << msg.getString("Mise à jour des valeurs") << " ...";
    
    for (size_t n=0; n<porosite.size(); ++n) {
        double cond = boost::numeric::ublas::inner_prod(lambda[n], phi);
        porosite[n] += cond;
    }
    if ( verbose ) std::cout << msg.getString("terminée") << '\n';
}

void Grille::lirePorosite(const string fichier)
{
    std::ifstream fin;
    
    fin.open( fichier.c_str() );
    if ( !fin.is_open() ) {
        std::cerr << msg.getString("Impossible d'ouvrir ") << fichier << std::endl;
        exit(1);
    }
    double dtmp;
    for (int n=0; n<9; ++n)
        fin >> dtmp; // nn.i nn.j nn.k   d.i d.j d.k   min.i min.j min.k
    
    std::vector<double> p;
    size_t I = 0;
    fin >> dtmp;  // on lit la 1re valeur
    while ( fin && I<porosite.size() ) {
        p.push_back( dtmp );
        fin >> dtmp;
        I++;
    }
    fin.close();
    if ( p.size() != porosite.size() ) {
        std::cerr << msg.getString("err2") << std::endl;
        exit(1);
    }
    else {
        for (size_t n=0; n<porosite.size(); ++n)
            porosite[n] = p[n];
    }
    simulerPorosite = false;
}
