#pragma once
#include "stdio.h"
#include <cmath>

template<typename T>
class Point{
	bool compare(T x, T y)
	{
		const double eps = 1E-10;
		if (abs(x - y) < eps)
			return true;
		else
			return false;
	}
public:
	T x;
	T y;
	T z;
	int index;

	Point()
	{this->x=0; this->y=0; this->z=0; this->index=0;}
	Point(double x, double y, double z)
	{this->x=x; this->y=y; this->z=z; this->index=0;}
	Point<T>& operator = ( const Point<T>& A )
	{
		this->x = A.x;
		this->y = A.y;
		this->z = A.z;
		this->index = A.index;
		return *this;
	}
	Point<T>& operator = (const T & A)
	{
		this->x = A;
		this->y = A;
		this->z = A;
		this->index = 0;
		return *this;
	}

	Point<T>& operator /= ( const T& A )
	{
		this->x = this->x / A;
		this->y = this->y / A;
		this->z = this->z / A;
		return *this;
	}
	Point<T>& operator /= (const Point<T>& A)
	{
		this->x = this->x / A.x;
		this->y = this->y / A.y;
		this->z = this->z / A.z;
		return *this;
	}
	Point<T>& operator *= ( const T& A )
	{
		this->x = this->x * A;
		this->y = this->y * A;
		this->z = this->z * A;
		return *this;
	}
	Point<T>& operator += ( const Point<T>& A )
	{
		this->x = this->x + A.x;
		this->y = this->y + A.y;
		this->z = this->z + A.z;
		return *this;
	}
	Point<T>& operator -= (const Point<T>& A)
	{
		this->x = this->x - A.x;
		this->y = this->y - A.y;
		this->z = this->z - A.z;
		return *this;
	}
	Point<T> operator + (const Point<T>& right)
	{
		return Point<T>(this->x + right.x, this->y + right.y, this->z + right.z);
	}
	Point<T> operator - (const Point<T>& right)
	{
		return Point<T>(this->x - right.x, this->y - right.y, this->z - right.z);
	}
	Point<T> operator / (const double& A)
	{
		return Point<T>(this->x/A, this->y/A, this->z/A);
	}
	Point<T> operator * (const double& A)
	{
		return Point<T>(this->x*A, this->y*A, this->z*A);
	}
	/*Point operator * (const std::vector<std::vector<double>> &A)
	{
		Point _X;
		_X.x = this->x*A[0][0] + this->y*A[0][1] + this->z*A[0][2];
		_X.y = this->x*A[1][0] + this->y*A[1][1] + this->z*A[1][2];
		_X.z = this->x*A[2][0] + this->y*A[2][1] + this->z*A[2][2];
		return _X;
	}
	Point operator * (const double A[3][3])
	{
		Point _X;
		_X.x = this->x*A[0][0] + this->y*A[0][1] + this->z*A[0][2];
		_X.y = this->x*A[1][0] + this->y*A[1][1] + this->z*A[1][2];
		_X.z = this->x*A[2][0] + this->y*A[2][1] + this->z*A[2][2];
		return _X;
	}*/

	//векторное произведение
	Point<T> operator ^ (const Point<T>& A)
	{
		return Point<T>(this->y*A.z-this->z*A.y, this->z*A.x-this->x*A.z, this->x*A.y-this->y*A.x);
	}
	//скалярное произведение
	double operator * (const Point<T>& A) 
	{
		return this->x*A.x+this->y*A.y+this->z*A.z;
	}
	
	bool operator == ( const Point<T> &A )
	{
		if( compare(this->x,A.x) && compare(this->y,A.y) && compare(this->z,A.z) )
			return true;
		return false;
	}
	bool operator != (const Point<T> &A)
	{
		if (!compare(this->x, A.x) || !compare(this->y, A.y) || !compare(this->z, A.z))
			return true;
		return false;
	}
	bool operator > ( const Point<T> &A )
	{
		if( compare(this->x,A.x) )
		{
			if( compare(this->y,A.y) )
			{
				if( this->z > A.z ) return true;
				if( this->z < A.z ) return false;
			}
			if( this->y > A.y ) return true;
			if( this->y < A.y ) return false;
		}
		if( this->x > A.x ) return true;
		if( this->x < A.x ) return false;
		return false;
	}
	bool operator < ( const Point<T> &A )
	{
		return !( *this > A );
	}
	bool operator >= ( const Point<T> &A )
	{
		return *this > A || *this == A;
	}
	bool operator <= ( const Point<T> &A )
	{
		return *this < A || *this == A;
	}
	Point<T> Pow( int k )
	{
		return Point<T>(pow(this->x,k), pow(this->y,k), pow(this->z,k));
	}
	bool include( Point<T> X )
	{
		double eps = 1E-3;
		int XX = int(abs(X.x-this->x) * 100000);
		int YY = int(abs(X.y-this->y) * 100000);
		int ZZ = int(abs(X.z-this->z) * 100000);
		//if ( abs(X.x-this->x)/abs(this->x)<=eps && abs(X.y-this->y)/abs(this->y)<=eps && abs(X.z-this->z)/abs(this->z)<=eps )
		if ( XX==0 && YY==0 && ZZ==0 )
			return true;
		return false;
	}

	void rotation(Point<T> center, double angleXY, double angleXZ, double angleYZ)
	{
		double _x, _y, _z;
		_x = this->x;
		_y = this->y;
		_z = this->z;

		_x = (this->x - center.x) * cos(angleXY) + (this->y - center.y) * sin(angleXY) + center.x;
		_y = -1 * (this->x - center.x) * sin(angleXY) + (this->y - center.y) * cos(angleXY) + center.y;
		this->x = _x;
		this->y = _y;
		this->z = _z;

		_x = (this->x - center.x) * cos(angleXZ) + (this->z - center.z) * sin(angleXZ) + center.x;
		_z = -1 * (this->x - center.x) * sin(angleXZ) + (this->z - center.z) * cos(angleXZ) + center.z;
		this->x = _x;
		this->y = _y;
		this->z = _z;

		_y = (this->y - center.y) * cos(angleYZ) + (this->z - center.z) * sin(angleYZ) + center.y;
		_z = -1 * (this->y - center.y) * sin(angleYZ) + (this->z - center.z) * cos(angleYZ) + center.z;
		this->x = _x;
		this->y = _y;
		this->z = _z;
	}

	void print()
	{
		printf_s("(%.2e, %.2e, %.2e)", this->x, this->y, this->z);
	}
};