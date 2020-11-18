// ------------------------------------------------------------------
// Clase Point - Punto en el plano
// ------------------------------------------------------------------

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cmath>
#include <vector>

#define PI 3.14159265
typedef float FLOAT;
//typedef double FLOAT;

using namespace std;

class Point {

public:
	/** coordinates */
	FLOAT x, y;

	Point (): x(0), y (0) {}
	Point(int ax, int ay) : x(ax), y(ay) {}
	Point(float ax, float ay) : x(ax), y(ay) {}
	Point (double ax, double ay): x (ax), y (ay) {}

/** Distance to other point */
	float dist(const Point& p) const
	{
		float dx = x - p.x;
		float dy = y - p.y;
		return sqrt (dx * dx + dy * dy);
	}
	
	bool operator== (const Point& p) const { return (x == p.x) && (y == p.y); }
	bool operator!= (const Point& p) const { return !(*this == p); }

	double distance(Point *otherPoint) const {
		return sqrt(pow((this->x - otherPoint->x), 2) + pow((this->y - otherPoint->y), 2));
	}

	/**
	* Calculate the distance between this point and the zero point
	* @return euclidean distance between this point and zero point
	*/
	double lengthFromZero() const {
		return sqrt(pow(this->x, 2) + pow(this->y, 2));
	}

	/**
	* Calculate the angle between two zero vectors
	* @param point2D other zero vector's point
	* @return angle between two vectors
	*/
	double phi(Point *point2D) {
		double vectorProduct = this->x * point2D->x + this->y * point2D->y;
		double lengthProduct = sqrt(this->x * this->x + this->y * this->y) *
			sqrt(point2D->x * point2D->x + point2D->y * point2D->y);
		return acos(vectorProduct / lengthProduct) * (180.0 / PI);
	}

};

inline ostream& operator<< (ostream& o, const Point& p) { return o << "(" << p.x << "," << p.y << ")"; }

inline istream& operator>> (istream& i, Point& p) { return i >> p.x >> p.y; }

#endif /* POINT_H */
