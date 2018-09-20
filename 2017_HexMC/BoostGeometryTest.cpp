#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>

namespace trans = boost::geometry::strategy::transform;
using boost::geometry::dsv;
using namespace boost::geometry;
BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian)
int main()
{


    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> point_type;

    point_type p1(1.0, 1.0);

    // Translate over (1.5, 1.5)
    point_type p2;
    trans::translate_transformer<double, 2, 2> translate(1.5, 1.5);
    boost::geometry::transform(p1, p2, translate);

    // Scale with factor 3.0
    point_type p3;
    trans::scale_transformer<double, 2, 2> scale(3.0);
    boost::geometry::transform(p1, p3, scale);

    // Rotate with respect to the origin (0,0) over 90 degrees (clockwise)
    point_type p4;
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(90.0);
    boost::geometry::transform(p1, p4, rotate);
	
	double box_line[][2] = {{2, 1}, {0, 1}, {0, -1}, {2, -1}, {2, 1}};
	model::polygon<model::d2::point_xy<double> > box;
	model::polygon<model::d2::point_xy<double> > box1;
	boost::geometry::append(box, box_line);
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate1(90.0);
    boost::geometry::transform(box, box1, rotate1);
	
	
    std::cout
        << "p1: " << dsv(p1) << std::endl
        << "p2: " << dsv(p2) << std::endl
        << "p3: " << dsv(p3) << std::endl
        << "p4: " << dsv(p4) << std::endl
		<< "box: " << dsv(box) << std::endl
		<< "box1: " << dsv(box1) << std::endl;

    return 0;
}