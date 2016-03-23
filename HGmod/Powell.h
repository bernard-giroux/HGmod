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

#ifndef __POWELL_H__
#define __POWELL_H__

#include <iostream>
#include <valarray>
#include <cmath>


template<typename T>
class Powell {
private:
	// routines adaptées du bouquin Numerical Recipes in C
	
	static const T gold;
	static const T glimit;
	static const T tiny;
	static const int itmax;
	static const T cgold;
	static const T zeps;
	static const T tol;
	
	static std::valarray<T> pcom;
	static std::valarray<T> xicom;
	static T (*nrfunc)(const std::valarray<T>&);
	
	static void mnbrak(T& ax, T& bx, T& cx, T& fa, T& fb, T& fc, T (*func)(T));
	static T brent(T ax, T bx, T cx, T (*f)(T), T tol, T& xmin);
	static T f1dim(T x);
	static void linmin(std::valarray<T>& p, std::valarray<T>& xi,
					   T& fret, T (*func)(const std::valarray<T>&));
	
public:
	static void powell(std::valarray<T>& p, std::valarray<T>& xi,
					   const T ftol, long& iter, T& fret,
					   T (*func)(const std::valarray<T>&));
	
    static T max(const T a, const T b);
	static T sign(const T a, const T b);
	static void shft(T& a, T& b, T& c, T& d);
	static T sqr(T a);
    

}; // end class

template <typename T>
class Brent {
    static const int itmax;
    static const T cgold;
	static const T zeps;

public:
    static T brent(T ax, T bx, T cx,
                   T (*f)(const T, const T, const T, const T, const T),
                   T tol, T& xmin,
                   std::valarray<T>& func_data);
};
#endif
