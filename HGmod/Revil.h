/*
 *  Revil.h
 *  HGmod
 *
 *  Created by Bernard Giroux on 08-01-21.
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

/*
*
*
 * @ARTICLE{revil97,
 * author = {Revil, A. and Glover, P. W. J.},
 * title = {Theory of ionic-surface electrical conduction in porous media },
 * journal = {Physical Review B},
 * year = {1997},
 * volume = {55},
 * pages = {1757--1773},
 * number = {3},
 * month = {Jan},
 * doi = {10.1103/PhysRevB.55.1757},
 * numpages = {16},
 * publisher = {American Physical Society}
 * }
 */


#ifndef __REVIL_H__
#define __REVIL_H__

template<typename T>
class Revil {
public:
    Revil(const Fluide<T> &fl, T F_, T f_, T L, T l)
    : fluide(fl), F(F_), f(f_), Lambda(L), lambda(l)
    {
        updateCtes();
    }
    
private:
	Fluide<T> fluide;
    // les 4 paramètres micro-géométriques
    T F;
    T f;
    T Lambda;
    T lambda;
    
    // les constantes de la fct G tilde (eq. 19)
    T a;
    T b;
    T c;
    T d;
    
    void updateCtes();
    T G(T X) { return (b + c*X + d*X*X)/(1.0 + a*X); }
};

template<typename T>
void Revil::updateCtes() {
    a = (2.0/(Lambda*F) - 1.0/f)/(0.5*lambda/f - 1.0/F);
    b = 1.0/F;
    c = (1.0 - lambda/Lambda)/(f - 0.5*lambda*F);
    d = a/f;
}

#endif