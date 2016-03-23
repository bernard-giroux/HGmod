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

#ifndef __STRUCTHGMOD_H__
#define __STRUCTHGMOD_H__

#include <valarray>
#include <vector>

#include "enum_bg.h"
#include "utils_bg.h"
#include "Geometrie.h"

struct paramCovariance
{
	double palier;
	double nbHurst;
	Utils_bg::typeCovariance typeCov;
	std::valarray<double> a;        // portées
	std::valarray<double> theta;    // angle anisotropie [degrés]
    
    paramCovariance() : palier(0.0001), nbHurst(0.0),
    typeCov(Utils_bg::EXPONENTIEL), a(std::valarray<double>(1.0,3)),
    theta(std::valarray<double>(0.0,3)) {}
};

struct paramConditionnement
{
    bool conditionne;
	std::vector<Utils_bg::double_ijk> coord;
	std::vector<double> data;
    
    paramConditionnement() : conditionne(false), coord(), data() {}
};

struct paramGenerateur
{
	int seed;                   // seed, nb aléatoire
	int L;                      // nb de ligne -> TBM
	Utils_bg::typeGenerateur typeGen;
	paramConditionnement pCond;
    
    paramGenerateur() : seed(1), L(1000), typeGen(Utils_bg::FFTMA), pCond() {}
};

template<typename T>
struct paramSaturation        // profil vertical de saturation
{
	T Smin;
	T Smax;
	T profNappe;
	T epaisseurTransition;
	T vG_alpha;
	T vG_m;
	T vG_n;
	bool vG;               // van Genuchten (80)  S(h) = [ 1 + (alpha h)^n ]^-m
    
    paramSaturation() : Smin(1.0), Smax(1.0), profNappe(3.0),
    epaisseurTransition(1.0), vG_alpha(0.0), vG_m(0.0), vG_n(0.0), vG(false) {}
};

template<typename TYPE>
struct paramPhysique
{
	
    bool KozenyCarman;              // calculer la perméabilité avec le modèle de Kozeny-Carman
    TYPE facteurPermea;             // facteur de forme pour le modèle de Kozeny-Carman si KozenyCarman==true
                                    // ou constante m de Pride si KozenyCarman==false
	TYPE phi_p;                     // porosité seuil de percolation - Kozeny-Carman
	TYPE TDS;                       // mg/l ou ppm
	TYPE temperature;               // Celsius
    
	TYPE mPorosite;                 // porosité totale moyenne
	TYPE Km;                        // const. diélectrique - matrice poreuse
	TYPE phi_s;                     // porosité intrinsèque du sable
	TYPE phi_c;                     // porosité intrinsèque de l'argile
	TYPE S_s;                       // surface spécifique intrinsèque du sable
	TYPE S_c;                       // surf. spécifique intrinsèque de l'argile
	TYPE T;                         // tortuosité
	TYPE sigma_c;                   // conductivité - argile
	TYPE Ks_c;                      // const. diélectrique statique - argile
	TYPE Ki_c;                      // const. diélectrique optique - argile
	TYPE fr_c;                      // fréq. relaxation - argile
	TYPE ccq_c;                     // Exposant Cole-Cole - argile
	TYPE m;                         // coefficient de cimentation
	TYPE n;                         // exposant de la Saturation
	TYPE permeabilite;              // relative
    TYPE V0;                        // P-wave velocity of rock matrix
    TYPE Vfl;                       // P-wave velocity of fluid
    
	std::vector<paramCovariance> pCov;    // covariance de la porosité
	std::vector<paramCovariance> sCov;    // covariance de la salinité
    paramConditionnement sData;           // données de salinité en mg/l ou ppm
	paramSaturation<TYPE> pSat;
    std::vector<double> granuloDiam;
    std::vector<double> granuloPassant;
    std::vector<double> granuloMasseVol;
    std::vector<double> granuloFactDepol;
    
    paramPhysique() : KozenyCarman(false), facteurPermea(8.0), phi_p(0.03),
    TDS(640.0), temperature(15.0), mPorosite(0.25), Km(5.0), phi_s(0.35),
    phi_c(0.45), S_s(2.65e5), S_c(5.0e8), T(2.5), sigma_c(-999999.0),
    Ks_c(30.0), Ki_c(20.0), fr_c(100.0e6), ccq_c(0.9), m(1.7), n(2.0),
    permeabilite(1.0), V0(6038.), Vfl(1500.), pCov(), sCov(), sData(), pSat(),
    granuloDiam(), granuloPassant(), granuloMasseVol(), granuloFactDepol()
    {}
};

struct paramSortie
{

	std::string fichierOut;            // nom du fichier contenant le modèle
	
	bool sauvePermittivite;
	bool sauveConductivite;
    bool sauvePermeaMag;
	bool sauvePermeabilite; // hydraulique
	bool sauvePorosite;
	bool sauveFractionArgile;
    bool sauveChargeabilite;
    bool sauveSalinite;
	bool sauveConductiviteComplexe;
    bool sauveVp;
    size_t  nf;  // nombre de freq du spectre
    
    paramSortie() : fichierOut(), sauvePermittivite(false),
    sauveConductivite(false), sauvePermeaMag(false), sauvePermeabilite(false),
	sauvePorosite(false), sauveFractionArgile(false), sauveChargeabilite(false),
    sauveSalinite(false), sauveConductiviteComplexe(false), sauveVp(false), nf(0) {}
    
};

struct paramGrille
{
	Utils_bg::double_ijk   d;          // incrément spatial selon x, y, z
	Utils_bg::double_ijk   min;        // coordonnées d'origine de la grille
	Utils_bg::double_ijk   max;        // coordonnées de la fin de la grille
	Utils_bg::size_t_ijk   nn;         // nombre de noeuds selon x, y, z
};

template<typename T>
struct paramEntree
{

	int nSimulations;
    bool simulerPorosite;
    std::string fichierPorosite;
	paramSortie pSortie;
	paramGrille  pGrille;
	paramPhysique<T> pPhys;
	paramGenerateur pGen;
	Geometrie* geom;
	T f;                         // fréquence de travail - radar [Hz]
    std::vector<T> frequences;   // fréquences du spectre de conductivité complexe [Hz]
    
    paramEntree() : nSimulations(1), simulerPorosite(true),
    fichierPorosite("porosity.dat"), pSortie(), pGrille(), pPhys(), pGen(),
    geom(0), f(100.0e6), frequences() {}
    
    paramEntree(const paramEntree<T>& p) : nSimulations(p.nSimulations),
    simulerPorosite(p.simulerPorosite), fichierPorosite(p.fichierPorosite),
    pSortie(p.pSortie), pGrille(p.pGrille), pPhys(p.pPhys), pGen(p.pGen),
    geom(p.geom), f(p.f), frequences(p.frequences) {}

    paramEntree<T>& operator=(const paramEntree<T>& p) {
        nSimulations = p.nSimulations;
        simulerPorosite = p.simulerPorosite;
        fichierPorosite = p.fichierPorosite;
        pSortie = p.pSortie;
        pGrille = p.pGrille;
        pPhys = p.pPhys;
        pGen = p.pGen;
        geom = p.geom;
        f = p.f;
        frequences = p.frequences;
    }

};

#endif
