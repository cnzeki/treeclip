#ifndef __TREE_CLIP_H__
#define __TREE_CLIP_H__
#include <vector>
#include "PolyStruct.h"

namespace  TreeClip {

	const float ALMOST_ZERO = 0.00001f;

	struct Point {
	public:
		float x;
		float y;

		Point(const Point &point) : x(point.x), y(point.y) {}
		Point(Point &point) : x(point.x), y(point.y) {}
		Point(float _x = 0., float _y = 0.) : x(_x), y(_y) {}

		//float x() const { return x; }

		//float y() const { return y; }

		void normalize()
		{
			float len = length();
			x /= len;
			y /= len;
		}

		Point normalized()
		{
			Point pt(x, y);
			pt.normalize();

			return pt;
		}

		void negate()
		{
			x = -x;
			y = -y;
		}

		float length() { return sqrt(x * x + y * y); }

		Point diff(Point pt) { return Point(x - pt.x, y - pt.y); }

		float distance(Point pt)
		{
			float dx = x - pt.x;
			float dy = y - pt.y;
			return sqrt(dx * dx + dy * dy);
		}

		float dot(Point pt) { return x*pt.x + y*pt.y; }

		Point orthogonal() { return Point(y, -x); }



		Point& operator-=(const Point &p) { x -= p.x; y -= p.y; return *this; }

		friend const Point operator+(const Point &p1, const Point &p2) { return Point(p1.x + p2.x, p1.y + p2.y); }
		friend const Point operator-(const Point &p1, const Point &p2) { return Point(p1.x - p2.x, p1.y - p2.y); }
		friend const Point operator*(const Point &p, float c) { return Point(p.x*c, p.y*c); }

		friend inline bool operator == (const Point &p1, const Point &p2);
		friend inline bool operator != (const Point &p1, const Point &p2);
	};

	inline bool operator == (const Point &p1, const Point &p2)
	{
		if (fabs(p1.x - p2.x) <= ALMOST_ZERO && fabs(p1.y - p2.y) <= ALMOST_ZERO)
		{
			return true;
		}

		return false;
	}

	inline bool operator != (const Point &p1, const Point &p2)
	{
		if (fabs(p1.x - p2.x) > ALMOST_ZERO || fabs(p1.y - p2.y) > ALMOST_ZERO)
		{
			return true;
		}

		return false;
	}

	typedef std::vector<Point>::iterator PointIterator;

	typedef std::vector< Point > Points;
	typedef std::vector< Point > Polygon;
	typedef std::vector< Polygon > Polygons;

	struct IntersectionPoint : public Point
	{
		bool entering;
		bool intersecting;
		bool used;

		//order of intersection in sequence
		int ord0, ord1;

		//index of start point of segment in each polygon
		int ind0, ind1;
		//distance of intersection to start point
		float dist0, dist1;
		//
		float step0, step1;
		//flag of intersection on segment,0:center,1:start point,2:end point
		int flag0, flag1;

		IntersectionPoint(float _x, float _y) : 
			Point(_x,_y), entering(false), intersecting(false), used(false), dist0(0.), dist1(0.), ord0(0), ord1(0), ind0(0), ind1(0), step0(0.), step1(0.), flag0(0), flag1(0) {}

		IntersectionPoint(float _x, float _y, bool _entering) :
			Point(_x, _y), entering(_entering), intersecting(true), used(false), dist0(0.), dist1(0.), ord0(0), ord1(0), ind0(0), ind1(0), step0(0.), step1(0.), flag0(0), flag1(0) {}

		IntersectionPoint(float _x, float _y, int _ind0, int _ind1) :
			Point(_x, _y), entering(false), intersecting(true), used(false), dist0(0.), dist1(0.), ord0(0), ord1(0), ind0(_ind0), ind1(_ind1), step0(_ind0), step1(_ind1), flag0(0), flag1(0) {}

		bool less0(IntersectionPoint &p2)
		{
			return step0 < p2.step0;
		}

		bool less1(IntersectionPoint &p2)
		{
			return step1 < p2.step1;
		}

		inline float calcDist(Point _startPoint, Point _endPoint)
		{
			float d0 = distance(_startPoint);
			float d1 = _endPoint.distance(_startPoint);

			return d0 / d1;
		}
	};

	typedef std::vector< IntersectionPoint > IntersectionPoints;

	class WeilerAtherton {
	public:
		WeilerAtherton();
		~WeilerAtherton();
		
		static IntersectionPoint doWalk(Polygon& walking_points, IntersectionPoints &inter_points, IntersectionPoint inter_point, bool entering, Polygon& pol);

		Polygons clipPoints(Polygon &subj_points, Polygon &clip_points, IntersectionPoints &subj_inter_points, IntersectionPoints &clip_inter_points);
	};
}

#endif /* __TREE_CLIP_H__ */
