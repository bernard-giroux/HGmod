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
#ifndef __GEOMETRIE_H__
#define __GEOMETRIE_H__

#include <vector>
#include <cmath>
#include <map>
#include <string>

enum typesGeometrie { CUBE, SPHERE, ELLIPSOIDE };

class TypesGeometrieNom
{
public:
    TypesGeometrieNom() : data()
   {
        data[CUBE] = "Cube";
        data[SPHERE] = "Sphere";
        data[ELLIPSOIDE] = "Ellipsoide";
   }
    std::string operator[](const typesGeometrie t) { return data[t]; }
    
private:
    std::map<typesGeometrie, std::string> data;
};

class Geometrie
{
public:
    Geometrie(typesGeometrie t) : numero_Medium(-1), tg(t) {}
	virtual ~Geometrie() {}
    void set_Numero_Medium(int n) { numero_Medium = n; }
    int get_Numero_Medium(){ return numero_Medium; }
    typesGeometrie get_typeGeometrie() const { return tg; }
    virtual bool a_l_interieur(double, double, double) const = 0;
    virtual double getXmin() const = 0;
    virtual double getXmax() const = 0;
    virtual double getYmin() const = 0;
    virtual double getYmax() const = 0;
    virtual double getZmin() const = 0;
    virtual double getZmax() const = 0;
private:
    int numero_Medium;
    typesGeometrie tg;
};

class Cube : public Geometrie
{
public:
    Cube(double xmi, double xma, double ymi, double yma,
         double zmi, double zma)
	: Geometrie(CUBE), xmin(xmi), xmax(xma), ymin(ymi), ymax(yma),
    zmin(zmi), zmax(zma) {}
    ~Cube() {}
    
    bool a_l_interieur(const double x, const double y, const double z) const
	{
		return ( x >= xmin && x <= xmax && y >= ymin && y <= ymax &&
				z >= zmin && z <= zmax );
	}
    
    double getXmin() const { return xmin; }
    double getXmax() const { return xmax; }
    double getYmin() const { return ymin; }
    double getYmax() const { return ymax; }
    double getZmin() const { return zmin; }
    double getZmax() const { return zmax; }
    
private:
    double xmin, xmax, ymin, ymax, zmin, zmax;
};

class Sphere : public Geometrie
{
public:
    Sphere(double x, double y, double z, double rayon)
	: Geometrie(SPHERE), x0(x), y0(y), z0(z), r(rayon) {}
    ~Sphere() {}
    
    bool a_l_interieur(const double x, const double y, const double z) const
	{
		return ( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) ) <= r;
	}
	
    
    double getXmin() const { return x0-r; }
    double getXmax() const { return x0+r; }
    double getYmin() const { return y0-r; }
    double getYmax() const { return y0+r; }
    double getZmin() const { return z0-r; }
    double getZmax() const { return z0+r; }
    
private:
    double x0, y0, z0, r;
};

class Ellipsoide : public Geometrie
{
public:
    Ellipsoide(double x, double y, double z, double r1, double r2, double r3,
               double alph, double thet)
	: Geometrie(ELLIPSOIDE), x0(x), y0(y), z0(z), a(r1), b(r2), c(r3),
    alpha(alph), theta(thet) {}
    ~Ellipsoide() {}
    
    bool a_l_interieur(const double x, const double y, const double z) const
	{
		// translation
		double xx = x-x0;
		double yy = y-y0;
		double zz = z-z0;
		// rotation
		double tmp = xx;
		// Dans le plan x-y, sens anti-horaire
		xx = xx*std::cos(theta)+yy*std::sin(theta);
		yy = -tmp*std::sin(theta)+yy*std::cos(theta);
		// Dans le plan x-z, sens anti-horaire
		tmp = xx;
		xx = xx*std::cos(alpha)+zz*std::sin(alpha);
		zz = -tmp*std::sin(alpha)+zz*std::cos(alpha);
		return ( (xx*xx)/(a*a) + (yy*yy)/(b*b) + (zz*zz)/(c*c) ) <= 1.;
	}
    
    double getXmin() const
	{
		double xx = x0-a;
		double yy = y0-b;
		double zz = z0-c;
		xx = xx*std::cos(theta)+yy*std::sin(theta);
		xx = xx*std::cos(alpha)+zz*std::sin(alpha);
		return xx;
	}
	
    double getXmax() const
	{
		double xx = x0+a;
		double yy = y0+b;
		double zz = z0+c;
		xx = xx*std::cos(theta)+yy*std::sin(theta);
		xx = xx*std::cos(alpha)+zz*std::sin(alpha);
		return xx;
	}
	
    double getYmin() const
	{
		double xx = x0-a;
		double yy = y0-b;
		//  double zz = z0-c;
		double tmp = xx;
		// Dans le plan x-y, sens anti-horaire
		//  xx = xx*std::cos(theta)+yy*std::sin(theta);
		yy = -tmp*std::sin(theta)+yy*std::cos(theta);
		return yy;
	}
	
    double getYmax() const
	{
		double xx = x0+a;
		double yy = y0+b;
		//  double zz = z0+c;
		double tmp = xx;
		// Dans le plan x-y, sens anti-horaire
		//  xx = xx*std::cos(theta)+yy*std::sin(theta);
		yy = -tmp*std::sin(theta)+yy*std::cos(theta);
		return yy;
	}
	
    double getZmin() const
	{
		double xx = x0-a;
		double yy = y0-b;
		double zz = z0-c;
		xx = xx*std::cos(theta)+yy*std::sin(theta);
		// Dans le plan x-z, sens anti-horaire
		double tmp = xx;
		xx = xx*std::cos(alpha)+zz*std::sin(alpha);
		zz = -tmp*std::sin(alpha)+zz*std::cos(alpha);
		return zz;
	}
	
    double getZmax() const
	{
		double xx = x0+a;
		double yy = y0+b;
		double zz = z0+c;
		xx = xx*std::cos(theta)+yy*std::sin(theta);
		// Dans le plan x-z, sens anti-horaire
		double tmp = xx;
		xx = xx*std::cos(alpha)+zz*std::sin(alpha);
		zz = -tmp*std::sin(alpha)+zz*std::cos(alpha);
		return zz;
	}
	
    
private:
    double x0, y0, z0, a, b, c, alpha, theta;
};

#endif
