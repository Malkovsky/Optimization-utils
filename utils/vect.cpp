#include <cstdio>
#include <cmath>
#include "vect.h"

	vect::vect(int n) {
		_size = n;
		#ifdef __RATIONAL_VECT_TYPE
		entry = new Rational[n];
		#endif
		#ifdef __DOUBLE_VECT_TYPE 
    	entry = new double[n];
		#endif
		for(int i = 0; i < n; ++i) {
			entry[i] = 0;
		}
	}

	vect_type &vect::operator[] (int i) {
		return entry[i];
	}
    
    vect_type vect::getEntry(int i) {
		if(i < 0 || i >= _size) {
			return -1;
		}
		return entry[i];
	}

	void vect::setEntry(int i, vect_type val) {
		if(i < 0 || i >= _size) {
			return;
		}
		entry[i] = val;
	}

	int vect::size() {
	    return _size;
	}
		
	vect_type vect::scalar_product(vect other) {
		if(_size != other._size) {
			return -1;
			//throw exception("Cannot multiply vectors, sizes are different!!!");
		}
		vect_type res = 0;    
		for(int i = 0; i < _size; ++i) {
			res += entry[i] * other.entry[i];
		}
		return res;
	}

	vect vect::operator + (vect other) {
    	if(_size != other._size) {
        	return *this;
    	}
    	vect* res = new vect(_size);
    	for(int i = 0; i < _size; ++i) {
    		res->entry[i] = entry[i] + other.entry[i];
   		}
   		return *res;		
	}

	vect vect::operator - (vect other) {
		if(_size != other._size) {
        	return *this;
    	}
    	vect* res = new vect(_size);
    	for(int i = 0; i < _size; ++i) {
    		res->entry[i] = entry[i] - other.entry[i];
   		}
   		return *res;		
	}

	vect vect::operator * (vect_type scalar) {
		vect* res = new vect(_size);
		for(int i = 0; i < _size; ++i) {
			res->entry[i] = scalar * entry[i];
		}
		return *res;
	}

	vect_type vect::operator * (vect other) {
		return scalar_product(other);
	}

	vect operator * (vect_type scalar, vect v) {
		return v * scalar;
	}

	vect vect::unities(int n) {
		vect res = vect(n);
		for(int i = 0; i < n; ++i) {
			res.entry[i] = 1;
		}
		return res;
	}

	vect vect::orth(int i, int n) {
		vect res = vect(n);
		res[i] = 1;
		return res;
	}

	vect_type vect::l22() {
		vect_type res = 0;
		for(int i = 0; i < _size; ++i) {
			res += entry[i] * entry[i];
		}
		return res; 
	}

	double vect::l2() {
		#ifdef __RATIONAL_VECT_TYPE
		return sqrt(l22().toDouble());
		#endif
		#ifdef __DOUBLE_VECT_TYPE
		return sqrt(l22());
		#endif
	}

	bool vect::operator == (vect other) {
		if(_size != other.size()) return false;
		bool ok = true;
		for(int  i = 0; i < _size; ++i ) {
			ok = ok && (entry[i] == other[i]);
		}
		return ok;
	}

	vect_type vect::abs_max_elem() {
		vect_type res = entry[0];
		for(int i = 0; i < _size; ++i) {
			if(res < fabs(entry[i])) res = fabs(entry[i]);			
		}
		return res;
	}

	void vect::print() {
		#ifdef __RATIONAL_VECT_TYPE
		for(int i = 0; i < _size; i++) 
			printf("%I64d/%I64d%c", entry[i].nom,  entry[i].denom, " \n"[i == _size - 1]);
		#endif
		#ifdef __DOUBLE_VECT_TYPE 
    	for(int i = 0; i < _size; i++) 
    		printf("%lf%c", entry[i], " \n"[i == _size - 1]);
    	#endif
    }

    void vect::err_print() {
		#ifdef __RATIONAL_VECT_TYPE
		for(int i = 0; i < _size; i++) 
			fprintf(stderr, "%I64d/%I64d%c", entry[i].nom,  entry[i].denom, " \n"[i == _size - 1]);
		#endif
		#ifdef __DOUBLE_VECT_TYPE 
    	for(int i = 0; i < _size; i++) 
    		fprintf(stderr, "%lf%c", entry[i], " \n"[i == _size - 1]);
    	#endif
    }

    std::string vect::tex_vert() {
    	std::string res;
    	res += "\\left(\\begin{array}{c}";
    	char buffer[62];
    	for(int i = 0; i < _size; ++i) {
    		//fprintf(stderr, "%I64d/%I64d\n", entry[i].nom, entry[i].denom);
    		#ifdef __RATIONAL_VECT_TYPE
    		//sprintf(buffer, "\\frac{%I64d}{%I64d}\\\\", entry[i].nom, entry[i].denom);
    		res += entry[i].tex(); 
    		res += "\\\\";
    		#endif
    		#ifdef __DOUBLE_VECT_TYPE
    		sprintf(buffer, "%lf\\\\", entry[i]);
    		res += buffer;
	    	#endif
	    }
	    res += "\\end{array}\\right)";
    	return res;
    }
    
	std::string vect::tex_gor(bool comma) {
    	std::string res;
    	res += "\\left(";
    	char buffer[52];
    	for(int i = 0; i < _size; ++i) {
    		#ifdef __RATIONAL_VECT_TYPE
    		//sprintf(buffer, "\\frac{%I64d}{%I64d}", entry[i].nom, entry[i].denom);
    		res += entry[i].tex(); 
    		#endif
    		#ifdef __DOUBLE_VECT_TYPE
    		sprintf(buffer, "%lf", entry[i]);
    		res += buffer;
    		#endif
    		if(comma) res += ",";
    		res += " \n"[i == _size - 1];
	    }
	    res += "\\right)";
	    return res;
    }
    
	