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


#ifndef __GRILLE_H__
#define __GRILLE_H__

#include <ostream>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "utils_bg_fct.h"
#include "structHGmod.h"
#include "Pride.h"
#include "Krigeage.h"
#include "Messages.h"
#include "ModelesMelange.h"

template<typename T, template<typename>class GEN, typename COV>
class MilieuPoreux;

extern bool verbose;
extern Messages msg;



class Grille
{
public:
	Grille(const paramGrille&, paramSortie&);
    
	template<typename T, template<typename>class GEN, typename COV>
        void remplire(const MilieuPoreux<T,GEN,COV>&,
                      const paramConditionnement&, const T,
                      const std::vector<T>&);
	template<typename T, template<typename>class GEN, typename COV>
		void remplireSalinite(const MilieuPoreux<T,GEN,COV>&);
	
    void sauveGrille(const std::string);
    void sauveSalinite(const std::string);
    template<typename T>
	void sauveData(const std::string&, const int, const std::vector<T>&);
    template<typename T>
	void sauveVTK(const std::string&, const int, const std::vector<T>&);
    void lirePorosite(const std::string);
    
    template<typename T, template<typename>class GEN, typename COV>
        void prepareCondPorosite(const MilieuPoreux<T,GEN,COV>&,
                                 const paramConditionnement&);
    
private:
        
    paramGrille p;
    paramSortie ps;
	// origine de la maille, avec pts pour conditions aux limites
	Utils_bg::double_ijk minm;
	// fin de la maille, avec pts pour conditions aux limites
	Utils_bg::double_ijk maxm;
    
    std::vector< boost::numeric::ublas::vector<double> > lambda;
	
	std::vector<double>                sigma;
	std::vector<std::complex<double> > permit;
	std::vector<double>                permea;
	std::vector<double>                porosite;
	std::vector<double>                cl;
	std::vector<double>                permeabilite;  // hydraulique
    std::vector<double>                chargeabilite;
    std::vector<double>                salinite;
    std::vector<double>                Vp;
    std::vector<short>                 facies;
	std::vector<std::vector<std::complex<double> > > sigmaCplx;
    bool saliniteVariable;
    bool simulerPorosite;
	
    template<typename T, template<typename>class GEN, typename COV>
		void remplirePorosite(const MilieuPoreux<T,GEN,COV>&);
	
    void conditionnePorosite(const paramConditionnement&);
};


template<typename T, template<typename>class GEN, typename COV>
void Grille::remplire(const MilieuPoreux<T,GEN,COV>& medium,
                      const paramConditionnement& pc, const T f,
                      const std::vector<T>& frequences)
{
    
    char barre[] = {'|', '/', '-', '\\'};
    if ( simulerPorosite ) {
        if ( pc.conditionne )
        {
            if ( ps.sauvePorosite==false ) {
                size_t nn = p.nn.i*p.nn.j*p.nn.k;
                porosite.resize(nn);
            }
            remplirePorosite(medium);
            conditionnePorosite(pc);
        }
    }
    double S, phi;
    
	if ( verbose ) std::cout << "\n" << msg.getString("Calcul des propriétés physiques ");
	if ( saliniteVariable ) {
		
		T temp = medium.getTemperature();
		T s = Pride<T>::ppm2mol_l(salinite[0]);
		T z[] = {1., -1.};
		T R[] = {116.e-12, 167.e-12}; // rayon ionique en m pour Na et Cl respectivement
		T eta = Utils_bg::viscositeEau(temp);  // viscosité en cP
                                               // mobilité des ions en N s/m (facteur 1000 pour passer de cP à kg/m.s)
		T b[] = {1000./(6.*Utils_bg::pi*eta*R[0]), 1000./(6.*Utils_bg::pi*eta*R[1])};
		T rho_w = Utils_bg::masseVolEau(salinite[0], temp);  // densité de l'eau
		T kappa = (ModelesMelange::MalmbergMaryott( temp ) +
                   ModelesMelange::corr_e_s(s))*Utils_bg::epsilon0;
        T pH = 7.0;
		Fluide<T> eau(s, 273.+temp, rho_w, eta, kappa, pH, z, b);
		
		for	( size_t ni=0, I=0; ni<p.nn.i; ++ni )
		{
			double x = minm.i + ni*p.d.i;
			if ( verbose ) {
				std::cout << "\r" << msg.getString("Calcul des propriétés physiques ") << barre[ni%4];
				std::cout.flush();
			}
			for ( size_t nj=0; nj<p.nn.j; ++nj )
			{
				double y = minm.j + nj*p.d.j;
				for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
				{
					double z = minm.k + nk*p.d.k;
					
					eau.setSalinite( Pride<T>::ppm2mol_l(salinite[I]) );
					if ( pc.conditionne || simulerPorosite==false )
						phi = porosite[I];
					else
						phi = medium.getPorositeTotale(x,y,z);
					
					if ( ps.sauvePorosite && pc.conditionne==false )
						porosite[I] = phi;
					
					if ( ps.sauveConductivite || ps.sauvePermittivite || 
                        ps.sauveConductiviteComplexe )
						S = medium.getSaturation(x, y, z);
					
					if (ps.sauveConductivite)	sigma[I] = medium.getSigma(S, phi, eau);
					if (ps.sauvePermittivite)   permit[I] = medium.getPermit(S, phi, f, eau);
					if (ps.sauveChargeabilite)  chargeabilite[I] = medium.getChargeabilite(S, phi);
                    
					if (ps.sauveFractionArgile) cl[I] = medium.getFractionArgile(phi);
					if (ps.sauvePermeabilite)   permeabilite[I] = medium.getPermeaHydro(phi);
					if (ps.sauvePermeaMag)      permea[I] = medium.getPermea(x, y, z);
                    if (ps.sauveConductiviteComplexe) {
                        for (size_t n=0; n<frequences.size(); ++n) {
                            sigmaCplx[I][n] = medium.getSigmaCplx(S, phi, frequences[n], eau);
                        }
                    }
                    if (ps.sauveVp) Vp[I] = medium.getVp(phi);
				}
			}
		}
	}
	else {
		for	( size_t ni=0, I=0; ni<p.nn.i; ++ni )
		{
			double x = minm.i + ni*p.d.i;
			if ( verbose ) {
				std::cout << "\r" << msg.getString("Calcul des propriétés physiques ")  << barre[ni%4];
				std::cout.flush();
			}
			for ( size_t nj=0; nj<p.nn.j; ++nj )
			{
				double y = minm.j + nj*p.d.j;
				for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
				{
					double z = minm.k + nk*p.d.k;
                    
					if ( pc.conditionne || simulerPorosite==false )
						phi = porosite[I];
					else
						phi = medium.getPorositeTotale(x,y,z);
                    
					if ( ps.sauvePorosite && pc.conditionne==false )
						porosite[I] = phi;
                    
					if ( ps.sauveConductivite || ps.sauvePermittivite || 
                        ps.sauveConductiviteComplexe )
						S = medium.getSaturation(x, y, z);
					if ( ps.sauveConductivite )
						sigma[I] = medium.getSigma(S, phi);
					if ( ps.sauvePermittivite )
						permit[I] = medium.getPermit(S, phi, f);
					if ( ps.sauveFractionArgile )
						cl[I] = medium.getFractionArgile(phi);
					if ( ps.sauvePermeabilite )
						permeabilite[I] = medium.getPermeaHydro(phi);
					if ( ps.sauvePermeaMag )
						permea[I] = medium.getPermea(x, y, z);
					if ( ps.sauveChargeabilite )
						chargeabilite[I] = medium.getChargeabilite(S, phi);
                    if ( ps.sauveConductiviteComplexe ) {
                        for (size_t n=0; n<frequences.size(); ++n) {
                            sigmaCplx[I][n] = medium.getSigmaCplx(S, phi, frequences[n]);
                        }
                    }
                    if (ps.sauveVp) Vp[I] = medium.getVp(phi);
				}
			}
		}
	}
    if ( verbose ) std::cout << " " << msg.getString("terminé") << "\n";
}

template<typename T, template<typename>class GEN, typename COV>
void Grille::remplireSalinite(const MilieuPoreux<T,GEN,COV>& medium)
{
	size_t nn = p.nn.i*p.nn.j*p.nn.k;
    if (salinite.size() != nn) salinite.resize(nn);
	std::valarray<double> x(3);
	
	char barre[] = {'|', '/', '-', '\\'};
	if ( verbose ) std::cout << "\n" << msg.getString("Calcul de la salinité ");
    for ( size_t ni=0, I=0; ni<p.nn.i; ++ni )
    {
        x[0] = minm.i + ni*p.d.i;
        if ( verbose ) {
            std::cout << "\r" << msg.getString("Calcul de la salinité ") << barre[ni%4];
            std::cout.flush();
        }
        for ( size_t nj=0; nj<p.nn.j; ++nj )
        {
            x[1] = minm.j + nj*p.d.j;
            for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
            {
                x[2] = minm.k + nk*p.d.k;
				salinite[I] = medium.getSalinite(x);
			}
		}
	}
	saliniteVariable = true;
}

template<typename T, template<typename>class GEN, typename COV>
void Grille::remplirePorosite(const MilieuPoreux<T,GEN,COV>& medium)
{
    char barre[] = {'|', '/', '-', '\\'};
    if ( verbose ) std::cout << "\n" << msg.getString("Calcul de la porosité ");
    for ( size_t ni=0, I=0; ni<p.nn.i; ++ni )
    {
        double x = minm.i + ni*p.d.i;
        if ( verbose ) {
            std::cout << "\r" << msg.getString("Calcul de la porosité ") << barre[ni%4];
            std::cout.flush();
        }
        for ( size_t nj=0; nj<p.nn.j; ++nj )
        {
            double y = minm.j + nj*p.d.j;
            for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
            {
                double z = minm.k + nk*p.d.k;
                porosite[I] = medium.getPorositeTotale(x,y,z);
            }
        }
    }
    if ( verbose ) std::cout << " " << msg.getString("terminé") << "\n";
}

template<typename T, template<typename>class GEN, typename COV>
void Grille::prepareCondPorosite(const MilieuPoreux<T,GEN,COV>&medium,
                                 const paramConditionnement& pc)
{
    if ( verbose )
        std::cout << "\n" << msg.getString("Conditionnement des porosités simulées (krigeage simple)");
    if ( verbose ) std::cout << "\n" << msg.getString("   Construction des matrices de krigeage");
    
    KrigeageSimple<COV> ks( medium.getGenerateur()->get_covariance_fct(),
                            &pc, 0.0);
    
    std::valarray<double> x1(3);
    
    char barre[] = {'|', '/', '-', '\\'};
    
    if ( verbose ) std::cout << "\n   " << msg.getString("Calcul des poids de krigeage ");
    for ( size_t ni=0, I=0; ni<p.nn.i; ++ni )
    {
        x1[0] = minm.i + ni*p.d.i;
        if ( verbose ) {
            std::cout << "\r   " << msg.getString("Calcul des poids de krigeage ") << barre[ni%4];
            std::cout.flush();
        }
        for ( size_t nj=0; nj<p.nn.j; ++nj )
        {
            x1[1] = minm.j + nj*p.d.j;
            for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
            {
                x1[2] = minm.k + nk*p.d.k;
                
                lambda.push_back( ks.getPoids(x1) );
            }
        }
    }
    if ( verbose ) std::cout << " " << msg.getString("terminé") << "\n";
}

template<typename T>
void Grille::sauveData(const std::string& fichier, const int ns,
                       const std::vector<T>& freq)
{
	if ( verbose ) {
		std::cout << "\nSauvegarde des données ... ";
		std::cout.flush();
	}
	char no[10];
	sprintf(no,"_%04d.dat",ns);
    std::ofstream fout((fichier+no).c_str());
    fout << '%';
    if (ps.sauvePorosite)       fout << '\t' << msg.getString("porosité");
	if (ps.sauvePermittivite)   fout << '\t'<<msg.getString("permittivité")<<"_Re\t"
		<<msg.getString("permittivité")<<"_Im";
    if (ps.sauveConductivite)   fout << '\t' << msg.getString("conductivité");
    if (ps.sauvePermeaMag)      fout << '\t' << msg.getString("perméabilité_magnétique");
    if (ps.sauvePermeabilite)   fout << '\t' << msg.getString("perméabilité_hydraulique");
    if (ps.sauveFractionArgile) fout << '\t' << msg.getString("fraction_argile");
    if (ps.sauveChargeabilite)  fout << '\t' << msg.getString("chargeabilité");
    if (ps.sauveConductiviteComplexe) {
        for (size_t n=0; n<freq.size(); ++n) {
            fout << '\t'<<msg.getString("conductivité")<<"_Re_"<< freq[n]<<"\t"
            <<msg.getString("conductivité")<<"_Im_"<< freq[n];
        }
    }
    if (ps.sauveVp)             fout << "\tVp";
    
    fout << '\n';
    
    for ( size_t ni=0, I=0; ni<p.nn.i; ++ni )
    {
        for ( size_t nj=0; nj<p.nn.j; ++nj )
        {
            for ( size_t nk=0; nk<p.nn.k; ++nk, I++ )
            {
                
                if (ps.sauvePorosite)       fout << '\t' << porosite[I];
                if (ps.sauvePermittivite)   fout << '\t' << permit[I].real() << '\t' << permit[I].imag();
                if (ps.sauveConductivite)   fout << '\t' << sigma[I];
                if (ps.sauvePermeaMag)      fout << '\t' << permea[I];
                if (ps.sauvePermeabilite)   fout << '\t' << permeabilite[I];
                if (ps.sauveFractionArgile) fout << '\t' << cl[I];
                if (ps.sauveChargeabilite)  fout << '\t' << chargeabilite[I];
                if (ps.sauveConductiviteComplexe) {
                    for (size_t n=0; n<freq.size(); ++n) {
                        fout << '\t'<< sigmaCplx[I][n].real() <<"\t"
                        << sigmaCplx[I][n].imag();
                    }
                }
                if (ps.sauveVp)             fout << '\t' << Vp[I];
                fout << '\n';
            }
        }
    }
    fout.close();
	if ( verbose ) std::cout << msg.getString("terminée") << '\n';
}

template<typename T>
void Grille::sauveVTK(const std::string& fichier, const int ns,
                       const std::vector<T>& freq)
{
#ifdef VTK
    
	if ( verbose ) {
		std::cout << "\nSauvegarde des données ... ";
		std::cout.flush();
	}
    
    vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
    for (size_t n=0; n<p.nn.i; ++n)
        xCoords->InsertNextValue(minm.i + n*p.d.i);
    vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
    for (size_t n=0; n<p.nn.j; ++n)
        yCoords->InsertNextValue(minm.j + n*p.d.j);
    vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
    for (size_t n=0; n<p.nn.k; ++n)
        zCoords->InsertNextValue(minm.k + n*p.d.k);
    
    vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
    rgrid->SetDimensions(p.nn.i, p.nn.j, p.nn.k);
    rgrid->SetXCoordinates(xCoords);
    rgrid->SetYCoordinates(yCoords);
    rgrid->SetZCoordinates(zCoords);
    
    int nc=0;
    if (ps.sauvePorosite) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("Porosity");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( porosite[ind] );
                }
            }
        }
        rgrid->GetPointData()->SetScalars( data );
        nc++;
    }
    if (ps.sauvePermittivite) {
        vtkSmartPointer<vtkDoubleArray> datar = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> datai = vtkSmartPointer<vtkDoubleArray>::New();
        datar->SetName("Permittivity - real part");
        datai->SetName("Permittivity - imaginary part");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    datar->InsertNextValue( permit[ind].real() );
                    datai->InsertNextValue( permit[ind].imag() );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( datar );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( datar );
        }
        rgrid->GetPointData()->AddArray( datai );
    }
    if (ps.sauveConductivite) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("Conductivity");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( sigma[ind] );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( data );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( data );
        }
    }
    if (ps.sauvePermeaMag) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("Magnetic permeability");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( permea[ind] );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( data );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( data );
        }
    }
    if (ps.sauvePermeabilite) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("Hydraulic permeability");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( permeabilite[ind] );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( data );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( data );
        }
    }
    if (ps.sauveFractionArgile) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("Clay fraction");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( cl[ind] );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( data );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( data );
        }
    }
    if (ps.sauveChargeabilite) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("Chargeability");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( chargeabilite[ind] );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( data );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( data );
        }
    }
    if (ps.sauveConductiviteComplexe) {
        for ( size_t nf=0; nf<freq.size(); ++nf ) {
            vtkSmartPointer<vtkDoubleArray> datar = vtkSmartPointer<vtkDoubleArray>::New();
            vtkSmartPointer<vtkDoubleArray> datai = vtkSmartPointer<vtkDoubleArray>::New();
            std::stringstream name;
            name << "Conductivity (" << freq[nf] << " Hz) - real part";
            datar->SetName( name.str().c_str() );
            name.seekp( 0, std::ios_base::beg );
            name << "Conductivity (" << freq[nf] << " Hz) - imaginary part";
            datai->SetName(name.str().c_str());
            for ( size_t nk=0; nk<p.nn.k; ++nk ) {
                for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                    for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                        size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                        datar->InsertNextValue( sigmaCplx[ind][nf].real() );
                        datai->InsertNextValue( sigmaCplx[ind][nf].imag() );
                    }
                }
            }
            if ( nc == 0 ) {
                rgrid->GetPointData()->SetScalars( datar );
                nc++;
            } else {
                rgrid->GetPointData()->AddArray( datar );
            }
            rgrid->GetPointData()->AddArray( datai );
        }
    }
    if (ps.sauveVp) {
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetName("P-wave velocity");
        for ( size_t nk=0; nk<p.nn.k; ++nk ) {
            for ( size_t nj=0; nj<p.nn.j; ++nj )  {
                for ( size_t ni=0; ni<p.nn.i; ++ni ) {
                    size_t ind = (ni*p.nn.j+nj)*p.nn.k+nk;
                    data->InsertNextValue( Vp[ind] );
                }
            }
        }
        if ( nc == 0 ) {
            rgrid->GetPointData()->SetScalars( data );
            nc++;
        } else {
            rgrid->GetPointData()->AddArray( data );
        }
    }

	char no[10];
	sprintf(no,"_%04d.vtr",ns);
    
    vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
    
    writer->SetFileName( (fichier+no).c_str() );
    writer->SetInputData( rgrid );
    writer->SetDataModeToBinary();
    writer->Update();

	if ( verbose ) std::cout << msg.getString("terminée") << '\n';
    
#endif
}



#endif

