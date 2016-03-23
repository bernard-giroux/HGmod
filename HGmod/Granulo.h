/*
 *  Granulo.h
 *  HGmod
 *
 *  Created by Bernard Giroux on 08-01-15.
 *  Copyright 2008 Bernard Giroux. All rights reserved.
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

#ifndef __GRANULO_H__
#define __GRANULO_H__

#include <cstdlib>
#include <iostream>
#include <vector>

class Granulo {
public:
    Granulo(const std::vector<double>& d, const std::vector<double>& p) :
    diametre(d), passant(p),
    masseVolumique(std::vector<double>(p.size(), 2.67)),
    fractionVolumique(std::vector<double>(p.size())),
    facteurDepolarisation(std::vector<double>(p.size(), 1./3.))
    {
        calculFractionVolumique();
    }

    Granulo(const std::vector<double>& d, const std::vector<double>& p,
            const std::vector<double>& m) :
    diametre(d), passant(p), masseVolumique(m),
    fractionVolumique(std::vector<double>(p.size())),
    facteurDepolarisation(std::vector<double>(p.size(), 1./3.))
    {
        calculFractionVolumique();
    }
    
    // retourne la fraction du volume des grains solides
    const std::vector<double>& getFractionVolumique() const {
        return fractionVolumique;
    }
    
    const std::vector<double>& getFacteurDepolarisation() const {
        return facteurDepolarisation;
    }
    
    void setFacteurDepolarisation( const std::vector<double>& f ) {
        if ( f.size() != facteurDepolarisation.size() ) {
            std::cerr << "Erreur : Granulo.h - taille incompatible" << std::endl;
            std::abort();
        }
        facteurDepolarisation = f;
    }
    
    void setMasseVolumique( const std::vector<double>& mv ) {
        if ( mv.size() != masseVolumique.size() ) {
            std::cerr << "Erreur : Granulo.h - taille incompatible" << std::endl;
            std::abort();
        }
        masseVolumique = mv;
        calculFractionVolumique();
    }
    
    size_t nData() const { return diametre.size(); }
    
    double getDiametre(size_t n) const { return diametre[n]; }
    
private:
    std::vector<double> diametre;  // doit être en ordre _décroissant_
    std::vector<double> passant;
    std::vector<double> masseVolumique;
    std::vector<double> fractionVolumique;
    std::vector<double> facteurDepolarisation;   // 1/3 -> sphere
    // 0 -> cylindre infiniment long (aiguille) d'axe parallèle au champ
    // 1 -> cylindre infiniment court (plaque) d'axe parallèle au champ
    
    void calculFractionVolumique();
    
}; // end class Granulo

#endif