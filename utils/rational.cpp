/**
 * @author Malkovsky Nikolay, april 2016
 *
 */

#include "rational.h"
#include <cstdlib>
#include <cstdio>

#ifdef __EUCLID_
long long gcd(long long a, long long b) {
    long long c;
	while (b != 0) {
		c = a % b;
		a = b;
		b = c;
	}
	if(a < 0) {
    	return -a;
    }
    return a;
}
#endif

#ifdef __TRIVIAL_GCD_
/**
 * Note that trivial gcd works only if a=b.
 */
long long gcd(long long a, long long b) {
	if (a == b || a == -b) return a > 0 ? a : -a;
	return 1;
}
#endif
    

Rational::Rational() {
	nom = 0;
	denom = 1;
}

Rational::Rational(long long a, long long b) {
	long long d = gcd(a, b);
	if(b == 0) {
		Rational();
   		return;
   	}
   	if(b < 0) {
   		b = -b;
   		a = -a;
   	}

   	nom = a / d;
   	denom = b / d;
}

Rational::Rational(long long a) {
	nom = a;
	denom = 1;
}

 Rational Rational::operator + ( Rational other) {
 	if(denom == other.denom) {
 		Rational ret;
 		ret.nom = nom + other.nom;
 		ret.denom = denom;
 		long long d = gcd(ret.nom, ret.denom);
 		ret.nom /= d;
 		ret.denom /= d;
 		return ret;
 	}
	long long c = nom * other.denom + other.nom * denom;
	long long f = denom * other.denom;
	long long d = gcd(c, f);
	Rational ret;
	ret.nom = c / d;
	ret.denom = f / d;
	return ret;
}

 Rational Rational::operator + ( int  other)  {
	Rational ret;
	ret.nom = nom + denom * other;
	ret.denom = denom;
	return ret;
}

 Rational Rational::operator + ( long long other)  {
	Rational ret;
	ret.nom = nom + denom * other;
	ret.denom = denom;
	return ret;
}

 Rational operator + ( long long other, Rational r)  {
	return r + other;
}

 Rational operator + ( int other, Rational r)  {
	return r + other;
}


Rational &Rational::operator += ( Rational other ) {
	return *this = *this + other;
}

Rational &Rational::operator += ( long long other ) {
	return *this = *this + other;
}

Rational &Rational::operator += ( int other ) {
	return *this = *this + other;
}

 Rational Rational:: operator - ( Rational other)  {
	if(denom == other.denom) {
 		Rational ret;
 		ret.nom = nom - other.nom;
 		ret.denom = denom;
 		long long d = gcd(ret.nom, ret.denom);
 		ret.nom /= d;
 		ret.denom /= d;
 		return ret;
 	}
	long long c = nom * other.denom - other.nom * denom;
	long long f = denom * other.denom;
	long long d = gcd(c, f);
	Rational ret;
	ret.nom = c / d;
	ret.denom = f / d;
	return ret;
}

 Rational Rational::operator - ( int other)  {
	Rational ret;
    ret.nom = nom - denom * other;
	ret.denom = denom;
	return ret;
}


 Rational Rational::operator - ( long long other)  {
	Rational ret;
	ret.nom = nom - denom * other;
	ret.denom = denom;
	return ret;
}

 Rational Rational:: operator - ()  {
	return Rational(-nom, denom);
}

 Rational operator - ( int other, Rational r)  {
	return (-r) + other;
}

 Rational operator - ( long long other, Rational r)  {
	return (-r) + other;
}

Rational &Rational::operator -= ( Rational other ) {
	return *this = *this - other;
}

Rational &Rational::operator -= ( long long other ) {
	return *this = *this - other;
}

Rational &Rational::operator -= ( int other ) {
	return *this = *this - other;
}


 Rational Rational:: operator * ( Rational other)  {
 	long long c = nom * (other.nom);
	//a = gcd(denom, other.denom);
	long long f = denom * (other.denom);
	long long d = gcd(c, f);
//	fprintf(stderr, "nom = %I64d denom = %I64d gcd = %I64d\n", c, f, d);
	Rational ret;
	ret.nom = c / d;
	ret.denom = f / d;
	return ret;
}

 Rational Rational:: operator * ( int other)  {
	long long d = gcd(denom, other);
	Rational ret;
	ret.nom = other / d * nom;
	ret.denom = denom / d;
	return ret;
}

 Rational Rational:: operator * ( long long other)  {
	long long d = gcd(denom, other);
	Rational ret;
	ret.nom = other / d * nom;
	ret.denom = denom / d;
	return ret;
}

 Rational  operator * ( long long other, Rational r)  {
	return r * other;
}

 Rational  operator * ( int other, Rational r)  {
	return r * other;
}

Rational &Rational::operator *= ( Rational other ) {
	return *this = *this * other;
}

Rational &Rational::operator *= ( long long other ) {
	return *this = *this * other;
}

Rational &Rational::operator *= ( int other ) {
	return *this = *this * other;
}


 Rational Rational::operator / ( Rational other)  {
	long long c = nom * other.denom;
	long long f = denom * other.nom;
	long long d = gcd(c, f);
	if(d < 0) d = -d;
	Rational ret;
//	fprintf(stderr, "nom = %I64d denom = %I64d gcd = %I64d\n", c, f, d);
	if(f > 0)  {
		ret.nom = c / d;
		ret.denom = f / d;
	} else {
		ret.nom = -c / d;
		ret.denom = -f / d;
	}
	return ret;
}

 Rational Rational::operator / ( int other)  {
	long long d = gcd(nom, other);
	Rational ret;
	if(d == 0) {
		fprintf(stderr, "d==0\n");
	}
	if(denom * other / d > 0) {
		ret.nom = nom / d;
		ret.denom = denom * other / d;
	} else {
		ret.nom = -nom / d;
		ret.denom = -denom * other / d;
	}
	return ret;
}

 Rational Rational::operator / ( long long other)  {
	long long d = gcd(nom, other);
	Rational ret;
	if(denom * other / d > 0) {
		ret.nom = nom / d;
		ret.denom = denom * other / d;
	} else {
		ret.nom = -nom / d;
		ret.denom = -denom * other / d;
	}
	return ret;
}

 Rational  operator / ( int other,  Rational r)  {
	long long d = gcd(r.nom, other);
	Rational ret;
	if(r.nom / d > 0) {
		ret.denom = r.nom / d;
		ret.nom = r.denom * other / d;
	} else {
		ret.denom = -r.nom / d;
		ret.nom = -r.denom * other / d;
	}
	return ret;
}

 Rational operator / ( long long other,  Rational r)  {
	long long d = gcd(r.nom, other);
	Rational ret;
	if(r.nom / d > 0) {
		ret.denom = r.nom / d;
		ret.nom = r.denom * other / d;
	} else {
		ret.denom = -r.nom / d;
		ret.nom = -r.denom * other / d;
	}
	return ret;
}

Rational &Rational::operator /= ( Rational other ) {
	return *this = *this / other;
}

Rational &Rational::operator /= ( long long other ) {
	return *this = *this / other;
}

Rational &Rational::operator /= ( int other ) {
	return *this = *this / other;
}

 bool Rational::operator < ( Rational other)  {
	return nom * other.denom < other.nom * denom;
}

 bool Rational::operator < ( long long other)  {
	return nom < other * denom;
}

 bool Rational::operator < ( int other)  {
	return nom < other * denom;
}

 bool operator < ( int other, Rational r)  {
	return other * r.denom < r.nom;
}

 bool operator <( long long other, Rational r)  {
	return other * r.denom < r.nom;
}

 bool Rational::operator > ( const Rational &other)  {
	return nom * other.denom > other.nom * denom;
}

 bool Rational::operator > ( const long long &other)  {
	return nom > other * denom;
}

 bool Rational::operator > ( const int &other)  {
	return nom > other * denom;
}

 bool operator > ( const int &other, const Rational &r)  {
	return other * r.denom > r.nom;
}

 bool operator >( const long long &other, const Rational &r)  {
	return other * r.denom > r.nom;
}

 bool Rational::operator <= ( const Rational &other)  {
	return nom * other.denom <= other.nom * denom;
}

 bool Rational::operator <= ( const long long &other)  {
	return nom <= other * denom;
}

 bool Rational::operator <= ( const int &other)  {
	return nom <= other * denom;
}

 bool operator <= ( const int &other, const Rational &r)  {
	return other * r.denom <= r.nom;
}

 bool operator <=( const long long &other, const Rational &r)  {
	return other * r.denom <= r.nom;
}

 bool Rational::operator >= ( const Rational &other)  {
	return nom * other.denom >= other.nom * denom;
}

 bool Rational::operator >= ( const long long &other)  {
	return nom >= other * denom;
}

 bool Rational::operator >= ( const int &other)  {
	return nom >= other * denom;
}

 bool operator >= ( const int &other, const Rational &r)  {
	return other * r.denom >= r.nom;
}

 bool operator >= ( const long long &other, const Rational &r)  {
	return other * r.denom >= r.nom;
}

 bool Rational::operator == ( const Rational &other)  {
	return nom * other.denom == other.nom * denom;
}

 bool Rational::operator == ( const long long &other)  {
	return nom == other * denom;
}

 bool Rational::operator == ( const int &other)  {
	return nom == other * denom;
}

 bool operator == ( const int &other, const Rational &r)  {
	return other * r.denom == r.nom;
}

 bool operator == ( const long long &other, const Rational &r)  {
	return other * r.denom == r.nom;
}


 bool Rational::operator != ( const Rational &other)  {
	return nom * other.denom != other.nom * denom;
}

 bool Rational::operator != ( const long long &other)  {
	return nom != other * denom;
}

 bool Rational::operator != ( const int &other)  {
	return nom != other * denom;
}

 bool operator != ( const int &other, const Rational &r)  {
	return other * r.denom != r.nom;
}

 bool operator != ( const long long &other, const Rational &r)  {
	return other * r.denom != r.nom;
}                           

double Rational::toDouble()  {
	return (double)nom / denom;
}

Rational fabs(Rational r) {
	return Rational(r.nom * (r.nom >= 0) - r.nom * (r.nom < 0), r.denom);
}

std::string Rational::tex() {
	std::string res;
	char buffer[20];
	if(denom == 1) {
		if(nom < 0) {
			res += "-";
			res += itoa(-nom, buffer, 10);
		} else {
			res += itoa(nom, buffer, 10);
		}
	} else {
		if(nom < 0)	{
			res += "-\\frac{";
			res += itoa(-nom, buffer, 10);
			res += "}{";
			res += itoa(denom, buffer, 10);
			res += "}";
		} else {
			res += "\\frac{";
			res += itoa(nom, buffer, 10);
			res += "}{";
			res += itoa(denom, buffer, 10);
			res += "}";
		}
	}
	return res;
}