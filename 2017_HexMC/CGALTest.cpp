#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <math.h>
#include <iostream>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Aff_transformation_2<K> Transformation; 
typedef CGAL::Vector_2<K> Vector; 
typedef CGAL::Direction_2<K> Direction; 


using std::cout; 
using std::endl;

int main1()
{
	
	
	
	
    Transformation t(CGAL::TRANSLATION, Vector(-2, 0));
    Transformation r(CGAL::ROTATION, sin(-M_PI/2), cos(-M_PI/2));

    Point points1[] = { Point(0, 0), Point(1, 1), Point(-1, 1), Point(-1, -1), Point(1, -1), Point(1, 1),};
    Point points2[] = { Point(0, 0), Point(1, 0), Point(0, 1)};
    Point points3[] = { Point(4, 0), Point(5, 0), Point(4, 1)};
    Polygon_2 pgn1(points1, points1+6);
    Polygon_2 pgn2(points2, points2+3);
    Polygon_2 pgn3(points3, points3+3);
    std::vector<Polygon_2> pgn_vec(0);
    pgn_vec.push_back(pgn1);
    pgn_vec.push_back(pgn2);
    pgn_vec.push_back(pgn3);
    bool flag = false;
    for(int i = 0; i < pgn_vec.size(); i++)
  	  if (pgn1 != pgn_vec[i])
  	    flag = flag || CGAL::do_intersect(pgn1, pgn_vec[i]);

    cout << "The polygon is " <<
      (flag ? "" : "not ") << "intersect." << endl;

 
    cout << "The polygon is " <<
      (CGAL::do_intersect(pgn1, pgn2) ? "" : "not ") << "intersect." << endl;
  
    cout << pgn_vec[0] << endl << CGAL::transform( t, pgn_vec[0]) << endl;
    cout << pgn_vec[0] << endl << CGAL::transform( r, pgn_vec[0]) << endl;
  
    return 0;
}