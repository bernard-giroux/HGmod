/*
 *  Messages.cpp
 *  HGmod
 *
 *  Created by Bernard Giroux on 06-11-17.
 *
 */
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

#include <cstdlib>
#include <iostream>
#include "Messages.h"

using namespace std;

Messages::Messages(string loc) : msg()
{
	if ( loc == "fr" ) constrFr();
	else if ( loc == "en" ) constrEn();
	else {
		cerr << "Messages: language not implemented" << endl;
		abort();
	}
}

void Messages::constrFr()
{
	msg["err1"] = "Erreur, générateur non implémenté.";
    msg["err2"] = "Format invalide - fichier de porosité";
    msg["err3"] = "Format invalide - il faut fournir un spectre pour modéliser la comductivité complexe";
	msg["Milieu Poreux"] = "Milieu Poreux";
	msg[" Simulation no "] = " Simulation no ";
	msg["Calcul des propriétés physiques "] = "Calcul des propriétés physiques ";
	msg["terminé"] = "terminé";
	msg["terminée"] = "terminée";
	msg["Calcul de la salinité "] = "Calcul de la salinité ";
	msg["Calcul de la porosité "] = "Calcul de la porosité ";
	msg["Conditionnement des porosités simulées (krigeage simple)"] = "Conditionnement des porosités simulées (krigeage simple)";
	msg["Construction des matrices de krigeage"] = "Construction des matrices de krigeage";
	msg["Calcul des poids de krigeage "] = "Calcul des poids de krigeage ";
	msg["Propriétés"] = "Propriétés";
	msg["Porosité moyenne"] = "Porosité moyenne";
	msg["Préparation de la grille"] = "Préparation de la grille";
	msg["La grille est de "] = "La grille est de ";
	msg["noeuds"] = "noeuds";
	msg["Pas: "] = "Pas: ";
	msg["Domaine modélisé"] = "Domaine modélisé";
	msg["Sauvegarde de la grille"] = "Sauvegarde de la grille";
	msg["grille"] = "grille";
	msg["Sauvegarde du modèle de salinité"] = "Sauvegarde du modèle de salinité";
	msg["salinité"] = "salinité";
	msg["porosité"] = "porosité";
	msg["permittivité"] = "permittivité";
	msg["conductivité"] = "conductivité";
	msg["perméabilité_magnétique"] = "perméabilité_magnétique";
	msg["perméabilité_hydraulique"] = "perméabilité_hydraulique";
	msg["fraction_argile"] = "fraction_argile";
	msg["chargeabilité"] = "chargeabilité";
	msg["Conditionnement des porosités simulées (krigeage simple)"] = "Conditionnement des porosités simulées (krigeage simple)";
	msg["Mise à jour des valeurs"] = "Mise à jour des valeurs";
	msg["Générateur: Classe de base"] = "Générateur: Classe de base";
	msg["Générateur: FFT-MA"] = "Générateur: FFT-MA";
	msg["Création du générateur FFT-MA"] = "Création du générateur FFT-MA";
	msg["Taille de la grille d'échantillonnage de la fct de covariance"] = "Taille de la grille d'échantillonnage de la fct de covariance";
	msg["Création du générateur terminée"] = "Création du générateur terminée";
	msg["Réinitialisation du générateur FFT-MA"] = "Réinitialisation du générateur FFT-MA";
	msg["Création du générateur par bandes tournantes"] = "Création du générateur par bandes tournantes";
	msg["Réinitialisation du générateur par bandes tournantes"] = "Réinitialisation du générateur par bandes tournantes";
	msg["Un générateur de modèles hydrogéophysiques"] = "Un générateur de modèles hydrogéophysiques";
	msg["Usage: "] = "Usage: ";
	msg["fichier.par"] = "fichier.par";
	msg["Affiche ce message"] = "Affiche ce message";
	msg["Spécification du fichier paramètres (obligatoire)"] = "Spécification du fichier paramètres (obligatoire)";
	msg["Affiche les messages informatifs"] = "Affiche les messages informatifs";
	msg["palier"] = "palier";
	msg["portées"] = "portées";
	msg["θ anisotropie"] = "θ anisotropie";
	msg["Modèle régional de covariance comportant"] = "Modèle régional de covariance comportant";
	msg["Création de la fonction de covariance effet pépite"] = "Création de la fonction de covariance effet pépite";
	msg["Covariance"] = "Covariance";
	msg["Effet pépite"] = "Effet pépite";
	msg["Création de la fonction de covariance de von Karman"] = "Création de la fonction de covariance de von Karman";
	msg["Modèle Von Karman"] = "Modèle Von Karman";
	msg["Nb de Hurst"] = "Nb de Hurst";
	msg["Création de la fonction de covariance exponentiel"] = "Création de la fonction de covariance exponentiel";
	msg["Modèle Exponentiel"] = "Modèle Exponentiel";
	msg["Création de la fonction de covariance sphérique"] = "Création de la fonction de covariance sphérique";
	msg["Modèle Sphérique"] = "Modèle Sphérique";
	msg["Création de la fonction de covariance gaussienne"] = "Création de la fonction de covariance gaussienne";
	msg["Modèle Gaussien"] = "Modèle Gaussien";
	msg["Création de la fonction de covariance hyperbolique"] = "Création de la fonction de covariance hyperbolique";
	msg["Modèle Hyperbolique"] = "Modèle Hyperbolique";
	msg["Création de la fonction de covariance Stable (α = 1/2)"] = "Création de la fonction de covariance Stable (α = 1/2)";
	msg["Modèle Stable (α = 1/2)"] = "Modèle Stable (α = 1/2)";
	msg["Erreur, modèle de covariance non implémenté"] = "Erreur, modèle de covariance non implémenté";
	msg["Impossible d'ouvrir "] = "Impossible d'ouvrir ";
	msg["Lecture du fichier paramètre: "] = "Lecture du fichier paramètre: ";
	msg["Température [°C]"] = "Température [°C]";
	msg["Permittivité relative du sable"] = "Permittivité relative du sable";
	msg["Porosité du sable pur"] = "Porosité du sable pur";
	msg["Porosité de l'argile pure"] = "Porosité de l'argile pure";
	msg["Surface spécifique du sable pur [m²/m³]"] = "Surface spécifique du sable pur [m²/m³]";
	msg["Surface spécifique de l'argile pure [m²/m³]"] = "Surface spécifique de l'argile pure [m²/m³]";
	msg["Tortuosité"] = "Tortuosité";
	msg["Modèle de salinité hétérogène"] = "Modèle de salinité hétérogène";
	msg["Salinité [ppm] (homogène)"] = "Salinité [ppm] (homogène)";
	msg["Permittivité relative statique de l'argile"] ="Permittivité relative statique de l'argile"; 
	msg["Permittivité relative optique de l'argile"] = "Permittivité relative optique de l'argile";
	msg["Fréquence de relaxation de l'argile [Hz]"] = "Fréquence de relaxation de l'argile [Hz]";
	msg["Exposant Cole-Cole - argile"] = "Exposant Cole-Cole - argile";
	msg["Fréquence d'évaluation de la permittivité [Hz]"] = "Fréquence d'évaluation de la permittivité [Hz]";
	msg["Modèle de saturation de van Genuchten appliqué"] = "Modèle de saturation de van Genuchten appliqué";
//	msg["Spécification du fichier de porosité (optionel)"] = "Spécification du fichier de porosité (optionel)";
    msg["Sauvegarde des fichiers en format VTK"] = "Sauvegarde des fichiers en format VTK";
}



void Messages::constrEn()
{
	msg["err1"] = "Error, generator not implemented.";
    msg["err2"] = "Invalid format - porosity file";
    msg["err3"] = "Invalid format - frequency spectrum needed to model complex conductivity";
	msg["Milieu Poreux"] = "Porous medium";
	msg[" Simulation no "] = " Simulation no ";
	msg["Calcul des propriétés physiques "] = "Computing physical properties ";
	msg["terminé"] = "done";
	msg["terminée"] = "done";
	msg["Calcul de la salinité "] = "Computing salinity ";
	msg["Calcul de la porosité "] = "Computing porosity ";
	msg["Conditionnement des porosités simulées (krigeage simple)"] = "Conditionning simulated porosities (simple kriging)";
	msg["Construction des matrices de krigeage"] = "Building kriging matrices";
	msg["Calcul des poids de krigeage "] = "Computing kriging weights ";
	msg["Propriétés"] = "Properties";
	msg["Porosité moyenne"] = "Mean porosity";
	msg["Préparation de la grille"] = "Preparing grid";
	msg["La grille est de "] = "The grid is ";
	msg["noeuds"] = "nodes";
	msg["Pas: "] = "Step size: ";
	msg["Domaine modélisé"] = "Domain";
	msg["Sauvegarde de la grille"] = "Saving grid";
	msg["grille"] = "grid";
	msg["Sauvegarde du modèle de salinité"] = "Saving salinity model";
	msg["salinité"] = "salinity";
	msg["porosité"] = "porosity";
	msg["permittivité"] = "permittivity";
	msg["conductivité"] = "conductivity";
	msg["perméabilité_magnétique"] = "magnetic_permeability";
	msg["perméabilité_hydraulique"] = "hydraulic_permeability";
	msg["fraction_argile"] = "clay_fraction";
	msg["chargeabilité"] = "chargeability";
	msg["Mise à jour des valeurs"] = "Updating values";
	msg["Générateur: Classe de base"] = "Generator: Base class";
	msg["Générateur: FFT-MA"] = "Generator: FFT-MA";
	msg["Création du générateur FFT-MA"] = "Creating FFT-MA generator";
	msg["Taille de la grille d'échantillonnage de la fct de covariance"] = "Size of covariance function sampling grid";
	msg["Création du générateur terminée"] = "Generator creation completed";
	msg["Réinitialisation du générateur FFT-MA"] = "Reinitializing FFT-MA generator";
	msg["Création du générateur par bandes tournantes"] = "Creating turning bands generator";
	msg["Réinitialisation du générateur par bandes tournantes"] = "Reinitializing turning bands generator";
	msg["Un générateur de modèles hydrogéophysiques"] = "A hydrogeophysical model builder";
	msg["Usage: "] = "Usage: ";
	msg["fichier.par"] = "file.par";
	msg["Affiche ce message"] = "Print this message";
	msg["Spécification du fichier paramètres (obligatoire)"] = "parameter file specification (mandatory)";
	msg["Affiche les messages informatifs"] = "Verbose";
	msg["palier"] = "sill";
	msg["portées"] = "ranges";
	msg["θ anisotropie"] = "θ anisotropy";
	msg["Modèle régional de covariance comportant"] = "Regional covariance model including";
	msg["Covariance"] = "Covariance";
	msg["Effet pépite"] = "Nugget Effect";
	msg["Modèle Von Karman"] = "Von Karman Model";
	msg["Nb de Hurst"] = "Hurst number";
	msg["Modèle Exponentiel"] = "Exponential Model";
	msg["Modèle Sphérique"] = "Sphérical Model";
	msg["Modèle Gaussien"] = "Gaussian Model";
	msg["Modèle Hyperbolique"] = "Hyperbolic Model";
	msg["Modèle Stable (α = 1/2)"] = "Stable Model (α = 1/2)";
	msg["Création de la fonction de covariance hyperbolique"] = "Creating hyperbolic covariance function";
	msg["Création de la fonction de covariance gaussienne"] = "Creating gaussian covariance function";
	msg["Création de la fonction de covariance effet pépite"] = "Creating nugget effect covariance function";
	msg["Création de la fonction de covariance de von Karman"] = "Creating von Karman covariance function";
	msg["Création de la fonction de covariance exponentiel"] = "Creating exponential covariance function";
	msg["Création de la fonction de covariance sphérique"] = "Creating sphérical covariance function";
	msg["Création de la fonction de covariance Stable (α = 1/2)"] = "Creating stable covariance function (α = 1/2)";
	msg["Erreur, modèle de covariance non implémenté"] = "Error, covariance model not implemented";
	msg["Impossible d'ouvrir "] = "Impossible to open ";
	msg["Lecture du fichier paramètre: "] = "Reading parameter file: ";
	msg["Température [°C]"] = "Temperature [°C]";
	msg["Permittivité relative du sable"] = "Relative permittivity of sand";
	msg["Porosité du sable pur"] = "Porosity of pur sand";
	msg["Porosité de l'argile pure"] = "Porosity of pur clay";
	msg["Surface spécifique du sable pur [m²/m³]"] = "Specific surface of pur sand [m²/m³]";
	msg["Surface spécifique de l'argile pure [m²/m³]"] = "Specific surface pur clay [m²/m³]";
	msg["Tortuosité"] = "Tortuosity";
	msg["Modèle de salinité hétérogène"] = "Heterogeneous salinity model";
	msg["Salinité [ppm] (homogène)"] = "Salinity [ppm] (homogeneous)";
	msg["Permittivité relative statique de l'argile"] ="Relative static permittivity of clay"; 
	msg["Permittivité relative optique de l'argile"] = "Relative optical permittivity of clay";
	msg["Fréquence de relaxation de l'argile [Hz]"] = "Relaxation frequency of clay [Hz]";
	msg["Exposant Cole-Cole - argile"] = "Cole-Cole exponent - clay";
	msg["Fréquence d'évaluation de la permittivité [Hz]"] = "Frequency of evaluation of permittivity [Hz]";
	msg["Modèle de saturation de van Genuchten appliqué"] = "Saturation model of van Genuchten applied";
    msg["Sauvegarde des fichiers en format VTK"] = "Save model files in VTK format";
}
