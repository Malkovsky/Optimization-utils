#ifndef __VECT
#define __VECT

//#define __DOUBLE_VECT_TYPE
#define __RATIONAL_VECT_TYPE
#include <string>


#ifdef __DOUBLE_VECT_TYPE
typedef double vect_type;
#endif

#ifdef __RATIONAL_VECT_TYPE
#include "rational.h"
typedef Rational vect_type;
#else
typedef double vect_type;
#endif

class vect {
	vect_type* entry;
	int _size;
public:
	vect(int n);
	vect_type &operator[] (int i);    
    vect_type getEntry(int i);
	void setEntry(int i, vect_type val);
	int size();		
	vect_type scalar_product(vect other);
	vect operator + (vect other);
	vect operator - (vect other);
	vect operator * (vect_type scalar);
	vect_type operator * (vect other);
	friend vect operator * (vect_type scalar, vect v);
	bool operator == (vect other);
	static vect unities(int n);
	static vect orth(int i, int n);
	vect_type l22();
	double l2();
	vect_type abs_max_elem();
	void print();
    void err_print();
    std::string tex_vert();
    std::string tex_gor(bool comma);
};

#endif