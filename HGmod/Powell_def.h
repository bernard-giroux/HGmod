/*
 *  Powell_def.h
 *  HGmod
 *
 *  Created by Bernard Giroux on 05-09-30.
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


#ifndef __POWELL_DEF_H__
#define __POWELL_DEF_H__

#include "Powell.h"


// Variables statiques

template<typename T>
const int Powell<T>::itmax = 5000;

template<typename T>
const T Powell<T>::gold   = 1.618034;

template<typename T>
const T Powell<T>::cgold  = 0.3819660;

template<typename T>
const T Powell<T>::glimit = 100.0;

template<typename T>
const T Powell<T>::tiny   = 1.0e-20;

template<typename T>
const T Powell<T>::zeps   = 1.0e-10;

template<typename T>
const T Powell<T>::tol    = 2.0e-4;

template<typename T>
std::valarray<T> Powell<T>::pcom;

template<typename T>
std::valarray<T> Powell<T>::xicom;

template<typename T>
T (*Powell<T>::nrfunc)(const std::valarray<T>&);


// Fonctions

template<typename T>
inline T Powell<T>::max(const T a, const T b)
{
	return a>b ? a : b;
}

template<typename T>
inline T Powell<T>::sign(const T a, const T b)
{
	return b>0.0 ? std::fabs(a) : -std::fabs(a);
}

template<typename T>
void Powell<T>::shft(T& a, T& b, T& c, T& d)
{
	a=b; b=c; c=d;
}

template<typename T>
inline T Powell<T>::sqr(T a)
{
	return a*a;
}


//
template<typename T>
void Powell<T>::mnbrak(T& ax, T& bx, T& cx, T& fa,
				   T& fb, T& fc, T (*func)(T))
{
	
	fa = (*func)(ax);
	fb = (*func)(bx);
	T dum;
	if (fb > fa) {
		shft(dum, ax, bx, dum);
		shft(dum, fb, fa, dum);
	}
	cx = bx + gold*(bx-ax);
	fc = (*func)(cx);
	while (fb > fc) {
		T fu;
		T r = (bx-ax)*(fb-fc);
		T q = (bx-cx)*(fb-fa);
		T u = bx-((bx-cx)*q - (bx-ax)*r)/(2.0*sign(max(std::fabs(q-r),tiny), q-r));
		T ulim = bx + glimit*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu=(*func)(u);
			if (fu < fc) {
				ax=bx;
				bx=u;
				fa=fb;
				fb=fu;
				return;
			} else if (fu > fb) {
				cx=u;
				fc=fu;
				return;
			}
			u = cx + gold*(cx-bx);
			fu = (*func)(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < fc) {
				T tmp = cx+gold*(cx-bx);
				shft(bx, cx, u, tmp);
				tmp = (*func)(u);
				shft(fb, fc, fu, tmp);
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u = ulim;
			fu = (*func)(u);
		} else {
			u = cx + gold*(cx-bx);
			fu = (*func)(u);
		}
		shft(ax, bx, cx, u);
		shft(fa, fb, fc, fu);
	}
}

//
template<typename T>
T Powell<T>::brent(T ax, T bx, T cx, T (*f)(T), T tol, T& xmin)
{
	T d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	T e=0.0;
	
	T a=((ax < cx) ? ax : cx);
	T b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (long iter=0; iter<itmax; iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*std::fabs(x)+zeps);
		if (std::fabs(x-xm) <= (tol2-0.5*(b-a))) {
			xmin=x;
			return fx;
		}
		if (std::fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=std::fabs(q);
			etemp=e;
			e=d;
			if (std::fabs(p) >= std::fabs(0.5*q*etemp) ||
				p <= q*(a-x) || p >= q*(b-x))
				d=cgold*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=sign(tol1,xm-x);
			}
		} else {
			d=cgold*(e=(x >= xm ? a-x : b-x));
		}
		u=(std::fabs(d) >= tol1 ? x+d : x+sign(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			shft(v,w,x,u);
			shft(fv,fw,fx,fu);
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	std::cerr << "Too many iterations in BRENT\n";
	xmin=x;
	return fx;
}


template<typename T>
T Powell<T>::f1dim(T x) {
	std::valarray<T> xt(pcom);
	for (size_t j=0; j<xicom.size(); j++) xt[j] += x*xicom[j];
	return (*nrfunc)(xt);
}

template<typename T>
void Powell<T>::linmin(std::valarray<T>& p, std::valarray<T>& xi,
				   T& fret, T (*func)(const std::valarray<T>&))
{
	pcom.resize(p.size());
	xicom.resize(xi.size());
	nrfunc=func;
	pcom = p;
	xicom = xi;
	T ax=0.0;
	T xx=1.0;
	T bx=2.0;
	T xmin,fx,fb,fa;
	mnbrak(ax, xx, bx, fa, fx, fb, f1dim);
	fret=brent(ax, xx, bx, f1dim, tol, xmin);
	for (size_t j=0; j<xi.size(); j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
}


// *** Modif par rapport a la routine powell de Numerical Recipes ***
// La routine powell prend une fonction avec un valarray de taille p.size()
// et minimise en changeant les sqrt(xi.size()) premiers elements de p

template<typename T>
void Powell<T>::powell(std::valarray<T>& p, std::valarray<T>& xi,
					   const T ftol, long& iter, T& fret,
					   T (*func)(const std::valarray<T>&))
{
	size_t n = static_cast<size_t>(sqrt(1.*xi.size()));
	
	std::valarray<T> pt = p;
	std::valarray<T> ptt(p.size());
	std::valarray<T> xit(n);
	
	fret = (*func)(p);
	
	for (iter=1; ; iter++) {
		T fp = fret;
		size_t ibig = 0;
		T del = 0.0;
		for (size_t i=0; i<n; i++) {
			for (size_t j=0; j<n; j++) xit[j]=xi[j*n+i];
			T fptt = fret;
			linmin(p,xit,fret,func);
			if (fabs(fptt-fret) > del) {
				del = fabs(fptt-fret);
				ibig = i;
			}
		}
		if (2.0*fabs(fp-fret) <= ftol*(fabs(fp)+fabs(fret))) {
			return;
		}
		if (iter == itmax) {
			std::cerr << "Too many iterations in routine POWELL\n";
			return;
		}
		for (size_t j=0; j<n; j++) {
			ptt[j] = 2.0*p[j] - pt[j];
			xit[j] = p[j] - pt[j];
			pt[j]  = p[j];
		}
		for (size_t j=n; j<p.size(); j++) {
			ptt[j] = p[j];
			pt = p[j];
		}
		T fptt=(*func)(ptt);
		if (fptt < fp) {
			T t=2.0*(fp-2.0*fret+fptt)*sqr(fp-fret-del)-del*sqr(fp-fptt);
			if (t < 0.0) {
				linmin(p, xit, fret, func);
				for (size_t j=0;j<n;j++) xi[j*n+ibig] = xit[j];
			}
		}
	}
}




template<typename T>
const int Brent<T>::itmax = 5000;

template<typename T>
const T Brent<T>::cgold  = 0.3819660;

template<typename T>
const T Brent<T>::zeps = std::numeric_limits<T>::epsilon()*1.e-3;


template<typename T>
T Brent<T>::brent(T ax, T bx, T cx,
                  T (*f)(const T, const T, const T, const T, const T),
                  T tol, T& xmin,
                  std::valarray<T>& func_data)
{
	T d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	T e=0.0;
	
	T a=((ax < cx) ? ax : cx);
	T b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x, func_data[0], func_data[1], func_data[2], func_data[3]);
	for (long iter=0; iter<itmax; iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*std::fabs(x)+zeps);
		if (std::fabs(x-xm) <= (tol2-0.5*(b-a))) {
			xmin=x;
			return fx;
		}
		if (std::fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=std::fabs(q);
			etemp=e;
			e=d;
			if (std::fabs(p) >= std::fabs(0.5*q*etemp) ||
				p <= q*(a-x) || p >= q*(b-x))
				d=cgold*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=Powell<T>::sign(tol1,xm-x);
			}
		} else {
			d=cgold*(e=(x >= xm ? a-x : b-x));
		}
		u=(std::fabs(d) >= tol1 ? x+d : x+Powell<T>::sign(tol1,d));
		fu=(*f)(u, func_data[0], func_data[1], func_data[2], func_data[3]);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			Powell<T>::shft(v,w,x,u);
			Powell<T>::shft(fv,fw,fx,fu);
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	std::cerr << "Too many iterations in BRENT\n";
	xmin=x;
	return fx;
}




#endif
