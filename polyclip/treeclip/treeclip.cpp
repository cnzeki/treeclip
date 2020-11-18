
#include <iostream>
#include "treeclip.h"

namespace TreeClip {
	WeilerAtherton::WeilerAtherton() {
	}

	WeilerAtherton::~WeilerAtherton() {
	}

	IntersectionPoint WeilerAtherton::doWalk(Polygon& walking_points, IntersectionPoints &inter_points, IntersectionPoint inter_point, bool entering, Polygon& pol)
	{
		int walking_size = walking_points.size();
		int inter_size = inter_points.size();

		auto it = std::find(inter_points.begin(), inter_points.end(), inter_point);
		int inter_start = it - inter_points.begin();
		int inter_end = (inter_start + 1) % inter_size;

		int poly_start, poly_end;
		int new_poly_start, new_poly_end;

		if (entering)
		{
			poly_start = inter_points[inter_start].ind0;
			poly_end = inter_points[inter_end].ind0;

			if (inter_point.flag0 == 0)
			{
				pol.emplace_back(inter_point.x, inter_point.y);
			}
		}
		else
		{
			poly_start = inter_points[inter_start].ind1;
			poly_end = inter_points[inter_end].ind1;
			if (inter_point.flag1 == 0)
			{
				pol.emplace_back(inter_point.x, inter_point.y);
			}
		}
		new_poly_start = (poly_start + 1) % walking_size;
		new_poly_end = (poly_end + 1) % walking_size;

		inter_points[inter_start].used = true;

		//if (inter_points[inter_end].entering)
		//{
		//	printf("enter:%.3f,%.3f\n", inter_points[inter_end].x, inter_points[inter_end].y);
		//}
		//else
		//{
		//	printf("exit:%.3f,%.3f\n", inter_points[inter_end].x, inter_points[inter_end].y);
		//}

		int i;
		if (new_poly_start < new_poly_end)
		{
			pol.insert(pol.end(), walking_points.begin() + new_poly_start, walking_points.begin() + new_poly_end);
		}
		else if(new_poly_start > new_poly_end)
		{
			pol.insert(pol.end(), walking_points.begin() + new_poly_start, walking_points.end());

			pol.insert(pol.end(), walking_points.begin(), walking_points.begin() + new_poly_end);
		}

		return inter_points[inter_end];
	}

	Polygons WeilerAtherton::clipPoints(Polygon &subj_points, Polygon &clip_points, IntersectionPoints &subj_inter_points, IntersectionPoints &clip_inter_points)
	{
		Polygons result;

		int i, j;
		IntersectionPoints::iterator it;
		for(i=0; i<subj_inter_points.size(); ++i)
		{
			Polygon pol;
			IntersectionPoint &start = subj_inter_points[i], next = start;
			if (start.used) continue;
			if(!start.entering) continue;

			//j = 0;
			do {
				//printf("----------------------i:%d,j:%d-----------------------\n", i, j);
				//printf("---subj walk start:%.3f,%.3f---\n", next.x, next.y);
				next = doWalk(subj_points, subj_inter_points, next, true, pol);

				//printf("---clip walk start:%.3f,%.3f---\n", next.x, next.y);
				next = doWalk(clip_points, clip_inter_points, next, false, pol);

				//j++;
			} while (next != start);

			result.push_back(pol);

			++i;
		}

		return result;
	}

}
