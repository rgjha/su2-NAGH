#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdlib.h>
#include "myrandom.h"

#include <ctime>
#include <time.h>

const int L = 4;
const int T = 1;
const int D = 2 ;

const int SITES = (L*T);

const int NCOLOR = 2;
const int NUMGEN = (NCOLOR*NCOLOR-1); 

// used for unitarity checks
const double GAUGETOL = 0.00000000000001;

const double PBC = 1.0;

const double SMALL = 0.00000000000001;
const int DEBUG_TIME = 1;

extern double BETA,KAPPA,C,DT,TIME;
extern int SWEEPS,GAP,START,THERM,READIN,SEED,TRAJECTORY_LENGTH;


class Complex{
	private:
		double re,im;
	public:
		Complex();
		Complex(double, double);
		double real(void) const;
		double imag(void) const;
		double norm(void);
		void print(void) const;
		friend ostream& operator<<(ostream&,Complex);
		friend istream& operator>>(istream&,Complex &);};
		
inline Complex conjug(const Complex &o1){return(Complex(o1.real(),-o1.imag()));}
inline Complex operator +(const Complex &o1, const Complex &o2){
        return(Complex(o1.real()+o2.real(),o1.imag()+o2.imag()));}
inline Complex operator -(const Complex &o1, const Complex &o2){
        return(Complex(o1.real()-o2.real(),o1.imag()-o2.imag()));}
inline Complex operator *(const Complex &o1, const Complex &o2){
        return(Complex(o1.real()*o2.real()-o1.imag()*o2.imag(),
                o1.real()*o2.imag()+o1.imag()*o2.real()));}
inline Complex operator *(const Complex &o1, const double o2){
        return(Complex(o1.real()*o2,o1.imag()*o2));}
inline Complex operator *(const double o1, const Complex &o2){
        return(Complex(o2.real()*o1,o2.imag()*o1));}

Complex operator /(const Complex &, const Complex &);
Complex pow(const Complex &, const int);


class Umatrix{
	private:
		Complex mat[NCOLOR][NCOLOR];
	public:
		Umatrix();
		Umatrix(int);
		Umatrix(Complex [NCOLOR][NCOLOR]);
		Complex get(const int,const int) const;
		void set(const int,const int,const Complex &);
		void print(void);
		friend ostream& operator<<(ostream &, Umatrix);
		friend istream& operator>>(istream &, Umatrix &);};
		
Umatrix operator +(const Umatrix &o1, const Umatrix &o2);
Umatrix operator -(const Umatrix &o1, const Umatrix &o2);
Umatrix operator *(const Umatrix &, const Umatrix &);
Umatrix operator *(const Umatrix &, const Complex &);
Umatrix operator *(const Complex &, const Umatrix &);
Umatrix operator *(const Umatrix &, const double);
Umatrix operator *(const double, const Umatrix &);
Umatrix exp(const Umatrix &u);
Umatrix Adj(const Umatrix &u);
Umatrix Cjg(const Umatrix &);
Complex Tr(const Umatrix &);
Umatrix gaussU(void);
Umatrix egaussU(void);
Umatrix cgaussU(void);

extern Umatrix Lambda[NUMGEN];

class Lattice_Vector{
private:
	int coords[D];
public:
	Lattice_Vector(void);
	Lattice_Vector(int);
	void set(int, int);
	int get(int) const;
	void print(void) const;
	};

Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x);
Lattice_Vector operator *(const int a, const Lattice_Vector &x);

double BC(const Lattice_Vector &x, const Lattice_Vector &y);
int loop_over_lattice(Lattice_Vector &, int &);

class Gauge_Field{
private:
	Umatrix link[SITES][D];
	
public:
    Gauge_Field(void);
	Gauge_Field(int);
	Umatrix get(const Lattice_Vector &, const int) const;
	void set(const Lattice_Vector &, const int, const Umatrix &);
	};
	
Gauge_Field Adj(const Gauge_Field &);
Gauge_Field Cjg(const Gauge_Field & );



#endif
