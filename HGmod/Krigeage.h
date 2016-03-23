/*
 *  Krigeage.h
 *  HGmod
 *
 *  Created by Bernard Giroux on 06-11-09.
 *  Copyright 2006 Bernard Giroux. All rights reserved.
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

#ifndef __KRIGEAGE_H__
#define __KRIGEAGE_H__

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


template<typename T>
class KrigeageOrdinaire
{
public:
    KrigeageOrdinaire(const T*, const paramConditionnement*);	
	~KrigeageOrdinaire() { delete Sigma; delete pivots; delete sigma; }
    
	double Z_est(const std::valarray<double>& x) const;
    boost::numeric::ublas::vector<double> getPoids(const std::valarray<double>& x) const;
private:
	size_t N;
	const paramConditionnement* data;
	const T* covar;

    boost::numeric::ublas::matrix<double> *Sigma;
    boost::numeric::ublas::permutation_matrix<double> *pivots;
    boost::numeric::ublas::vector<double> *sigma;
};

template<typename T>
KrigeageOrdinaire<T>::KrigeageOrdinaire(const T* covar_,
                                        const paramConditionnement* data_)
: covar(covar_), data(data_)
{
	N = data->coord.size();
	
	Sigma = new boost::numeric::ublas::matrix<double>(N+1, N+1);
	std::valarray<double> x1(3), x2(3);
	
	for (size_t n1=0; n1<N; n1++)
	{
		x1[0] = data->coord[n1].i;
		x1[1] = data->coord[n1].j;
		x1[2] = data->coord[n1].k;
		for (size_t n2=n1; n2<N; n2++)
		{
			x2[0] = data->coord[n2].i;
			x2[1] = data->coord[n2].j;
			x2[2] = data->coord[n2].k;
			(*Sigma)(n1,n2) = covar->C(x1, x2);
			if (n1!=n2) (*Sigma)(n2,n1) = (*Sigma)(n1,n2);
		}
		(*Sigma)(n1,N) = 1.0;
		(*Sigma)(N,n1) = 1.0;
	}
	(*Sigma)(N,N) = 0.0;
	
	pivots = new boost::numeric::ublas::permutation_matrix<double>(N+1);
	sigma = new boost::numeric::ublas::vector<double>(N+1);
	boost::numeric::ublas::lu_factorize(*Sigma, *pivots);
	(*sigma)[N] = 1.0;
}

template<typename T>
double KrigeageOrdinaire<T>::Z_est(const std::valarray<double>& x) const
{
	std::valarray<double> x2(3);
	for ( size_t n=0; n<N; ++n ) {
		x2[0] = data->coord[n].i;
		x2[1] = data->coord[n].j;
		x2[2] = data->coord[n].k;
		(*sigma)[n] = covar->C(x, x2);
	}
	boost::numeric::ublas::vector<double> lambda( *sigma );
    
	boost::numeric::ublas::lu_substitute(*Sigma, *pivots, lambda);
	double z = 0.0;
	for (size_t n=0; n<N; ++n) z += lambda[n]*data->data[n];
	return z;
}


template<typename T>
boost::numeric::ublas::vector<double> KrigeageOrdinaire<T>::getPoids(const std::valarray<double>& x) const
{
	std::valarray<double> x2(3);
	for ( size_t n=0; n<N; ++n ) {
		x2[0] = data->coord[n].i;
		x2[1] = data->coord[n].j;
		x2[2] = data->coord[n].k;
		(*sigma)[n] = covar->C(x, x2);
	}
	boost::numeric::ublas::vector<double> lambda( *sigma );
	boost::numeric::ublas::lu_substitute(*Sigma, *pivots, lambda);
	return lambda;
}



template<typename T>
class KrigeageSimple
{
public:
    KrigeageSimple(const T*, const paramConditionnement*, const double);
	~KrigeageSimple() { delete Sigma; delete pivots; delete sigma; delete Z; }
    
	double Z_est(const std::valarray<double>& x) const;
    boost::numeric::ublas::vector<double> getPoids(const std::valarray<double>& x) const;
private:
    size_t N;
	const paramConditionnement* data;
	const T* covar;
    boost::numeric::ublas::matrix<double> *Sigma;
    boost::numeric::ublas::permutation_matrix<double> *pivots;
    boost::numeric::ublas::vector<double> *sigma;
    boost::numeric::ublas::vector<double> *Z;
};

template<typename T>
KrigeageSimple<T>::KrigeageSimple(const T* covar_,
                                  const paramConditionnement* data_,
                                  const double moyenne)
: covar(covar_), data(data_)
{
	N = data->coord.size();
	
	Sigma = new boost::numeric::ublas::matrix<double>(N, N);
	std::valarray<double> x1(3), x2(3);
	
	for (size_t n1=0; n1<N; n1++)
	{
		x1[0] = data->coord[n1].i;
		x1[1] = data->coord[n1].j;
		x1[2] = data->coord[n1].k;
		for (size_t n2=n1; n2<N; n2++)
		{
			x2[0] = data->coord[n2].i;
			x2[1] = data->coord[n2].j;
			x2[2] = data->coord[n2].k;
			(*Sigma)(n1,n2) = covar->C(x1, x2);
			if (n1!=n2) (*Sigma)(n2,n1) = (*Sigma)(n1,n2);
		}
	}
	
	pivots = new boost::numeric::ublas::permutation_matrix<double>(N);
	sigma = new boost::numeric::ublas::vector<double>(N);
	boost::numeric::ublas::lu_factorize(*Sigma, *pivots);
    
    Z = new boost::numeric::ublas::vector<double>(N);
    for (size_t n=0; n<N; ++n) (*Z)[n] = data->data[n]-moyenne;
}

template<typename T>
double KrigeageSimple<T>::Z_est(const std::valarray<double>& x) const
{
	std::valarray<double> x2(3);
	for ( size_t n=0; n<N; ++n ) {
		x2[0] = data->coord[n].i;
		x2[1] = data->coord[n].j;
		x2[2] = data->coord[n].k;
		(*sigma)[n] = covar->C(x, x2);
	}
    boost::numeric::ublas::vector<double> lambda( *sigma );
	boost::numeric::ublas::lu_substitute(*Sigma, *pivots, lambda);
	return boost::numeric::ublas::inner_prod(lambda, *Z);
}


template<typename T>
boost::numeric::ublas::vector<double> KrigeageSimple<T>::getPoids(const std::valarray<double>& x) const
{
	std::valarray<double> x2(3);
	for ( size_t n=0; n<N; ++n ) {
		x2[0] = data->coord[n].i;
		x2[1] = data->coord[n].j;
		x2[2] = data->coord[n].k;
		(*sigma)[n] = covar->C(x, x2);
	}
    boost::numeric::ublas::vector<double> lambda( *sigma );
	boost::numeric::ublas::lu_substitute(*Sigma, *pivots, lambda);
	return lambda;
}


#endif
