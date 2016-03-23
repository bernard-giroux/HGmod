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

#ifndef __HGMOD_H__
#define __HGMOD_H__

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>

#include "structHGmod.h"
#include "Messages.h"

extern bool verbose;
extern Messages msg;

namespace HGmod {
	
	void print_usage (std::ostream& stream, int exit_code);
	
	std::string traite_arg(int argc, char *argv[]);
	
	template<typename T>
		void setDefaut(paramEntree<T>& p)	{
			p.nSimulations = 1;
            p.simulerPorosite = true;
            p.fichierPorosite = "porosity.dat";
			p.pGrille.min.i = 0.0;
			p.pGrille.min.j = 0.0;
			p.pGrille.min.k = 0.0;
			p.pGrille.nn.i = 101;
			p.pGrille.nn.j = 101;
			p.pGrille.nn.k = 101;
			p.pGrille.d.i = 0.1;
			p.pGrille.d.j = 0.1;
			p.pGrille.d.k = 0.1;
			p.pGrille.max.i = nan("0");
			p.pGrille.max.j = nan("0");
			p.pGrille.max.k = nan("0");
            
			p.pPhys.KozenyCarman = false;
			p.pPhys.facteurPermea = 8.0;
			p.pPhys.phi_p = 0.03;
			p.pPhys.mPorosite = 0.25;
			p.pPhys.TDS = 640.0;
			p.pPhys.temperature = 15.0;
			p.pPhys.Km = 5.0;
			p.pPhys.phi_s = 0.35;
			p.pPhys.phi_c = 0.45;
			p.pPhys.S_s = 2.65e5;
			p.pPhys.S_c = 5.0e8;
			p.pPhys.T = 2.5;
			p.pPhys.sigma_c = -999999.0;
			p.pPhys.Ks_c = 30.0;
			p.pPhys.Ki_c = 20.0;
			p.pPhys.fr_c = 100.e6;
			p.pPhys.ccq_c = 0.9;
			p.f = 100.e6;
			p.pPhys.m = 1.7;
            p.pPhys.n = 2.0;
            
            p.pPhys.pSat.Smin = 1.0;
            p.pPhys.pSat.Smax = 1.0;
            p.pPhys.pSat.profNappe = 3.0;
            p.pPhys.pSat.epaisseurTransition = 1.0;
			p.pPhys.pSat.vG = false;
			p.pPhys.pSat.vG_alpha = 0.0;
			p.pPhys.pSat.vG_m = 0.0;
			p.pPhys.pSat.vG_n = 0.0;
            
            p.pPhys.permeabilite = 1.0;
			
			p.pGen.typeGen = Utils_bg::FFTMA;
			p.pGen.pCond.conditionne = false;
			p.pGen.L = 1000;
            
			p.pPhys.pCov.resize(1);
			p.pPhys.pCov[0].typeCov = Utils_bg::EXPONENTIEL;
			p.pPhys.pCov[0].palier = 0.0001;
			p.pPhys.pCov[0].a.resize(3);
			p.pPhys.pCov[0].a = 1.0;
			p.pPhys.pCov[0].theta.resize(3);
			p.pPhys.pCov[0].theta = 0.0;
            
            p.pPhys.sData.conditionne = false;
            
            p.frequences = std::vector<T>();
            
		}
	
	template<typename T>
		void cleanup(paramEntree<T>& p) {
			delete p.geom;
		}
	
	
	template<typename T>
		void get_param(const std::string fichier, paramEntree<T>& p) {
            setDefaut(p);
			std::ifstream fin;
			
			fin.open( fichier.c_str() );
			if ( !fin.is_open() ) {
				std::cerr << msg.getString("Impossible d'ouvrir ") << fichier << std::endl;
				exit(1);
			}
			
			if ( verbose )
				std::cout << "\n" << msg.getString("Lecture du fichier paramètre: ") << fichier << std::endl;
			
			char valeur[100];
			char parametre[200];
			std::string par;
			std::istringstream sin( valeur );
            
			while (!fin.eof()) {
				fin.get(valeur, 100, '#');
				if (strlen(valeur) == 0) {
					fin.clear();
					fin.getline(parametre, 200);
					continue;
				}
				fin.get(parametre, 200, ',');
				par = parametre;
				
				if (par.find("Nombre de simulations") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.nSimulations;
				}
                
				if (par.find("Fichier Sortie") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pSortie.fichierOut;
				}
                
                if (par.find("Simuler Porosite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.simulerPorosite;
                }
                
				if (par.find("Fichier Porosite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.fichierPorosite;
				}
                
                if (par.find("Sauvegarde Permittivite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pSortie.sauvePermittivite;
                }
                
                if (par.find("Sauvegarde ConductiviteElectrique") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pSortie.sauveConductivite;
                }
                
                if (par.find("Sauvegarde PermeabiliteHydraulique") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pSortie.sauvePermeabilite; // hydraulique
                }
                
                if (par.find("Sauvegarde PermeabiliteMagnetique") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pSortie.sauvePermeaMag;
                }
                
                if (par.find("Sauvegarde Porosite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pSortie.sauvePorosite;
				}
                
                if (par.find("Sauvegarde FractionArgile") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pSortie.sauveFractionArgile;
                }
                
                if (par.find("Sauvegarde Chargeabilite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pSortie.sauveChargeabilite;
                }
                
                if (par.find("Sauvegarde Salinite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pSortie.sauveSalinite;
                }
                
                if (par.find("Sauvegarde ConductiviteComplexe") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pSortie.sauveConductiviteComplexe;
                }
                
                if (par.find("Sauvegarde Vp") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pSortie.sauveVp;
                }
				//
				// Grille
				//
				else if (par.find("Grille dx") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pGrille.d.i;
					sin >> p.pGrille.d.j;
					sin >> p.pGrille.d.k;
				}
				else if (par.find("Grille nNoeuds") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pGrille.nn.i;
					sin >> p.pGrille.nn.j;
					sin >> p.pGrille.nn.k;
				}
				else if (par.find("Grille x min") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pGrille.min.i;
					sin >> p.pGrille.min.j;
					sin >> p.pGrille.min.k;
				}
				else if (par.find("Grille x max") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pGrille.max.i;
					sin >> p.pGrille.max.j;
					sin >> p.pGrille.max.k;
				}
				//
				// Propriétés physiques
				//
                else if (par.find("Permeabilite Kozeny") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.KozenyCarman;
				}
                else if (par.find("Permeabilite facteur de forme") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.facteurPermea;
				}
                else if (par.find("Permeabilite porosite percolation") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.phi_p;
				}
				else if (par.find("Porosite moyenne") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.mPorosite;
				}
				else if (par.find("TDS") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.TDS;
				}
				else if (par.find("Temperature") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.temperature;
				}
				else if (par.find("Constante dielectrique matrice") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.Km;
				}
				
				else if (par.find("Porosite sable pur") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.phi_s;
				}
				else if (par.find("Porosite argile pure") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.phi_c;
				}
				else if (par.find("Surface specifique sable pur") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.S_s;
				}
				else if (par.find("Surface specifique argile pure") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.S_c;
				}
				else if (par.find("Tortuosite") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.T;
				}
				
				else if (par.find("conductivite argile") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.sigma_c;
				}
				else if (par.find("dielectrique statique argile") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.Ks_c;
				}
				else if (par.find("dielectrique optique argile") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.Ki_c;
				}
				else if (par.find("frequence relaxation argile") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.fr_c;
				}
				else if (par.find("exposant Cole-Cole argile") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.ccq_c;
				}
				
				else if (par.find("Saturation exposant") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.n;
				}
				else if (par.find("Saturation min") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.Smin;
				}				
				else if (par.find("Saturation max") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.Smax;
				}
				else if (par.find("Saturation profondeur nappe") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.profNappe;
				}
				else if (par.find("Saturation epaisseur frange capillaire") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.epaisseurTransition;
				}
				else if (par.find("Saturation van Genuchten") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.vG;
				}
				else if (par.find("Saturation vG alpha") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.vG_alpha;
				}
				else if (par.find("Saturation vG m") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.vG_m;
				}
				else if (par.find("Saturation vG n") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.pSat.vG_n;
				}
				else if (par.find("P-wave velocity - rock matrix") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.V0;
				}
				else if (par.find("P-wave velocity - fluid") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.Vfl;
				}
                
				else if (par.find("Covariance porosite - modele") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.pCov.size() < noCov ) p.pPhys.pCov.resize(noCov);
					std::string tmp;  sin >> tmp;
					if ( tmp.find("exp") < 20 )
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::EXPONENTIEL;
					else if ( tmp.find("Karman") < 20 ) 
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::VON_KARMAN;
					else if ( tmp.find("sph") < 20 )
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::SPHERIQUE;
					else if ( tmp.find("gauss") < 20 )
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::GAUSSIEN;
					else if ( tmp.find("pepite") < 20 )
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::PEPITE;
					else if ( tmp.find("hyperbol") < 20 )
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::HYPERBOLIQUE;
					else if ( tmp.find("stable") < 20 )
						p.pPhys.pCov[noCov-1].typeCov = Utils_bg::STABLE;
				}
				else if (par.find("Covariance porosite - palier") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.pCov.size() < noCov ) p.pPhys.pCov.resize(noCov);
					sin >> p.pPhys.pCov[noCov-1].palier;
				}
				else if (par.find("Covariance porosite - portees") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.pCov.size() < noCov ) p.pPhys.pCov.resize(noCov);
                    p.pPhys.pCov[noCov-1].a.resize(3);
					for (int n=0; n<3; ++n) sin >> p.pPhys.pCov[noCov-1].a[n];
				}
				else if (par.find("Covariance porosite - angles anisotropie") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.pCov.size() < noCov ) p.pPhys.pCov.resize(noCov);
                    p.pPhys.pCov[noCov-1].theta.resize(3);
					for (int n=0; n<3; ++n) sin >> p.pPhys.pCov[noCov-1].theta[n];
				}
				else if (par.find("Covariance porosite - Nombre de Hurst") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.pCov.size() < noCov ) p.pPhys.pCov.resize(noCov);
					sin >> p.pPhys.pCov[noCov-1].nbHurst;
				}
				
                else if (par.find("Covariance salinite - modele") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.sCov.size() < noCov ) p.pPhys.sCov.resize(noCov);
					std::string tmp;  sin >> tmp;
					if ( tmp.find("exp") < 20 )
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::EXPONENTIEL;
					else if ( tmp.find("Karman") < 20 ) 
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::VON_KARMAN;
					else if ( tmp.find("sph") < 20 )
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::SPHERIQUE;
					else if ( tmp.find("gauss") < 20 )
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::GAUSSIEN;
					else if ( tmp.find("pepite") < 20 )
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::PEPITE;
					else if ( tmp.find("hyperbol") < 20 )
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::HYPERBOLIQUE;
					else if ( tmp.find("stable") < 20 )
						p.pPhys.sCov[noCov-1].typeCov = Utils_bg::STABLE;
				}
				else if (par.find("Covariance salinite - palier") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.sCov.size() < noCov ) p.pPhys.sCov.resize(noCov);
					sin >> p.pPhys.sCov[noCov-1].palier;
				}
				else if (par.find("Covariance salinite - portees") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.sCov.size() < noCov ) p.pPhys.sCov.resize(noCov);
                    p.pPhys.sCov[noCov-1].a.resize(3);
					for (int n=0; n<3; ++n) sin >> p.pPhys.sCov[noCov-1].a[n];
				}
				else if (par.find("Covariance salinite - angles anisotropie") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.sCov.size() < noCov ) p.pPhys.sCov.resize(noCov);
                    p.pPhys.sCov[noCov-1].theta.resize(3);
					for (int n=0; n<3; ++n) sin >> p.pPhys.sCov[noCov-1].theta[n];
				}
				else if (par.find("Covariance salinite - Nombre de Hurst") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					size_t noCov;
					sin >> noCov;
					if ( p.pPhys.sCov.size() < noCov ) p.pPhys.sCov.resize(noCov);
					sin >> p.pPhys.sCov[noCov-1].nbHurst;
				}
                else if (par.find("Salinite - krigeage") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pPhys.sData.conditionne;
                }
                else if (par.find("Salinite - Point de mesure") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					Utils_bg::double_ijk tmp;
					sin >> tmp.i;
					sin >> tmp.j;
					sin >> tmp.k;
					p.pPhys.sData.coord.push_back(tmp);
					double tt;
					sin >> tt;
					p.pPhys.sData.data.push_back(tt);
				}
                
				else if (par.find("Porosite - generateur") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					std::string tmp;  sin >> tmp;
					if ( tmp.find("tbm") < 20 || tmp.find("TBM") < 20 )
						p.pGen.typeGen = Utils_bg::TBM;
					else if ( tmp.find("fftma") < 20 || tmp.find("FFTMA") < 20 ) 
						p.pGen.typeGen = Utils_bg::FFTMA;
				}
				else if (par.find("Porosite - seed") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pGen.seed;
				}
				else if (par.find("Covariance porosite - nombre de lignes TBM") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pGen.L;
				}
                else if (par.find("Porosite - conditionnement") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    sin >> p.pGen.pCond.conditionne;
                }
				else if (par.find("Porosite - Point de conditionnement") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					Utils_bg::double_ijk tmp;
					sin >> tmp.i;
					sin >> tmp.j;
					sin >> tmp.k;
					p.pGen.pCond.coord.push_back(tmp);
					double tt;
					sin >> tt;
					p.pGen.pCond.data.push_back(tt);
				}
				
				else if (par.find("Permeabilite relative") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.pPhys.permeabilite;
				}
                
                else if (par.find("Granulometrie - diametre des grains") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    double tmp;
                    while (sin >> tmp ) {
                        p.pPhys.granuloDiam.push_back(tmp);
                    }
                }
                else if (par.find("Granulometrie - fraction passante") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    double tmp;
                    while (sin >> tmp ) {
                        p.pPhys.granuloPassant.push_back(tmp);
                    }
                }
                else if (par.find("Granulometrie - masse volumique") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    double tmp;
                    while (sin >> tmp ) {
                        p.pPhys.granuloMasseVol.push_back(tmp);
                    }
                }
                else if (par.find("Granulometrie - facteur de depolarisation") < 200) {
                    sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
                    double tmp;
                    while (sin >> tmp ) {
                        p.pPhys.granuloFactDepol.push_back(tmp);
                    }
                }
                
                else if (par.find("Frequence de travail") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					sin >> p.f;
				}
                else if (par.find("Spectre de frequences") < 200) {
					sin.str( valeur ); sin.seekg(0, std::ios_base::beg); sin.clear();
					double tmp;
                    while (sin >> tmp ) {
                        p.frequences.push_back(tmp);
                    }
				}
				
				fin.getline(parametre, 200);
			}
			fin.close();
            
            if ( p.pSortie.sauveConductiviteComplexe && p.frequences.size()==0 ) {
                std::cerr << msg.getString("err3") << std::endl;
                exit(1);
            }
            p.pSortie.nf = p.frequences.size();
			
            if ( p.simulerPorosite ) {
                if (finite(p.pGrille.max.i)==0) {
                    p.pGrille.max.i = p.pGrille.min.i + (p.pGrille.nn.i-1)*p.pGrille.d.i;
                    p.pGrille.max.j = p.pGrille.min.j + (p.pGrille.nn.j-1)*p.pGrille.d.j;
                    p.pGrille.max.k = p.pGrille.min.k + (p.pGrille.nn.k-1)*p.pGrille.d.k;
                }
                else {
                    p.pGrille.nn.i = (unsigned int) ( 1.00000001+(p.pGrille.max.i-p.pGrille.min.i) / p.pGrille.d.i );
                    p.pGrille.nn.j = (unsigned int) ( 1.00000001+(p.pGrille.max.j-p.pGrille.min.j) / p.pGrille.d.j );
                    p.pGrille.nn.k = (unsigned int) ( 1.00000001+(p.pGrille.max.k-p.pGrille.min.k) / p.pGrille.d.k );
                }
            }
            else {
                
                fin.open( p.fichierPorosite.c_str() );
                if ( !fin.is_open() ) {
                    std::cerr << msg.getString("Impossible d'ouvrir ") << p.fichierPorosite << std::endl;
                    exit(1);
                }
                
                fin >> p.pGrille.nn.i;
                fin >> p.pGrille.nn.j;
                fin >> p.pGrille.nn.k;
                
                fin >> p.pGrille.d.i;
                fin >> p.pGrille.d.j;
                fin >> p.pGrille.d.k;
                
                fin >> p.pGrille.min.i;
                fin >> p.pGrille.min.j;
                fin >> p.pGrille.min.k;
                
                p.pGrille.max.i = p.pGrille.min.i + (p.pGrille.nn.i-1)*p.pGrille.d.i;
                p.pGrille.max.j = p.pGrille.min.j + (p.pGrille.nn.j-1)*p.pGrille.d.j;
                p.pGrille.max.k = p.pGrille.min.k + (p.pGrille.nn.k-1)*p.pGrille.d.k;
                
                fin.close();
            }
            
            p.geom = new Cube(p.pGrille.min.i, p.pGrille.max.i,
                              p.pGrille.min.j, p.pGrille.max.j,
                              p.pGrille.min.k, p.pGrille.max.k);
		}
};

#endif
