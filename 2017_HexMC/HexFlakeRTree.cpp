#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/foreach.hpp>

using namespace boost::geometry;
using std::cout; 
using std::endl;
BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian);

typedef model::d2::point_xy<double> Point_2;
typedef model::polygon<Point_2> Polygon_2;
typedef model::box<Point_2> Box_2;
typedef std::pair<Box_2, Polygon_2> Value;
typedef index::rtree< Value, index::quadratic<1000> > Lattice;


// Класс для генератора случайных чисел
class RandomGenerator{

private:
	unsigned long x, y, z, w;

private:
	unsigned long xor128()
	{
		unsigned long t;
		t = (x^(x<<11));
		x = y;
		y = z;
		z = w;
		return(w = (w^(w>>19))^(t^(t>>8)));
	};
	
// Конструктор. Инициализирует генератор.
public:
	RandomGenerator()
	{
		x = (unsigned long)(time(NULL));
		y = 362436069;
		z = 521288629;
		w = 88675123;
	};
// Следующее псевдослучайное число
public:
	double NextDouble()
	{
		return ((double)xor128()) / ULONG_MAX;
	};
};

// Одна итерация построения гексафлека
std::vector<Point_2> HexFlakeStep(std::vector<Point_2> &old_line)
{	
    std::vector<Point_2> new_line(0);
	size_t old_line_len = old_line.size();
	new_line.push_back(old_line[0]);
    for(int i = 0; i < old_line_len - 1; i++ )
    {
		Point_2 p;
		Point_2 p1;
		Point_2 p2;
		Point_2 p3;
		Point_2 p4;
		Point_2 p5;
        p1 = old_line[i];
        p2 = old_line[i+1];
        double dpx = (get<0>(p2) - get<0>(p1))/3.0;
        double dpy = (get<1>(p2) - get<1>(p1))/3.0;
		new_line.push_back(p1);
        p3 = Point_2(get<0>(p1) + dpx, get<1>(p1) + dpy);
        p4 = Point_2(get<0>(p1) + 2*dpx, get<1>(p1) + 2*dpy);
        double alpha = M_PI/3.0;
        p = Point_2((get<0>(p2) - get<0>(p1))/3.0, (get<1>(p2) - get<1>(p1))/3.0);
        p5 = Point_2(get<0>(p)*cos(alpha) - get<1>(p)*sin(alpha), get<1>(p)*cos(alpha) + get<0>(p)*sin(alpha));
        p5 = Point_2(get<0>(p5) + get<0>(p3), get<1>(p5) + get<1>(p3));
        new_line.push_back(p3);
        new_line.push_back(p5);
        new_line.push_back(p4);
	}
	new_line.push_back(old_line[old_line_len-1]);
    return  new_line;
};


// Функция сохранения флеков решетки lattice в файл "hexaflake.txt"

/*
void SaveLattice(Lattice &lattice, std::string file_name)
{

	std::ofstream out_file;
	out_file.open (file_name + ".txt");
	std::string str = "";
	BOOST_FOREACH(Value const& v, lattice)
		str += wkt<Polygon_2>(v.second);
	out_file << str;
	out_file.close();
}
*/

// Двигаем полигон с индексом pgn_index
Polygon_2 MovePolygon(Polygon_2 &pgn, double max_translation, double max_angle, RandomGenerator &rnd, 
					  double xmin, double xmax, double ymin, double ymax)
{
	Polygon_2 ret, ret1;
	Point_2 center;
	centroid(pgn, center);
	if( rnd.NextDouble() > 0.5)
	{
		double dx = -get<0>(center);
		double dy = -get<1>(center);
		double angle = max_angle*(1 - 2*rnd.NextDouble());
	    strategy::transform::translate_transformer<double, 2, 2> t1(dx, dy);
		transform(pgn, ret, t1);
	    strategy::transform::rotate_transformer<boost::geometry::radian, double, 2, 2> r(angle);
		transform(ret, ret1, r);
	    strategy::transform::translate_transformer<double, 2, 2> t2(-dx, -dy);
		transform(ret1, ret, t2);
	}
	else{
		double x = get<0>(center);
		double y = get<1>(center);
		double dx = max_translation*(1 - 2*rnd.NextDouble());
		double dy = max_translation*(1 - 2*rnd.NextDouble());
		if (x + dx > xmax)
			dx = dx - (xmax - xmin);
		else if (x + dx < xmin)
			dx = dx + (xmax - xmin);
		if (y + dy > ymax)
			dy = dy - (ymax - ymin);
		else if (y + dy < ymin)
			dy = dy + (ymax - ymin);
	    strategy::transform::translate_transformer<double, 2, 2> t(dx, dy);
		transform(pgn, ret, t);
	}
	return ret;
}

// Проверяем, что полигон pgn_probe, полученный из полигона с индексом pgn_index, никого не перекрывает
bool CheckPolygonIntersect(Lattice &lattice, Polygon_2 &pgn_probe)
{

	bool flag = false;
	Box_2 box;
	envelope(pgn_probe, box);
	std::vector<Value> result_s;
	lattice.query(index::intersects(box), std::back_inserter(result_s));
	//lattice.query(index::nearest(box, 10), std::back_inserter(result_s));
	BOOST_FOREACH(Value const& v, result_s)
		if(overlaps(v.second, pgn_probe))
		{
			flag = true;
			break;
		}
	return flag;
}

void Run(Lattice &lattice, 
		int steps_count, int nx, int ny, double xmin, double xmax, double ymin, double ymax)
{
	RandomGenerator rnd;
	Polygon_2 pgn_probe, pgn;
	// Максимальные смещения
	double max_translation = 0.1;
	double max_angle = 0.1;
	double accept = 0;
	for(int i = 0; i < steps_count; i++)
	{
		cout << "Time step " << i << endl;
		auto bt = std::chrono::high_resolution_clock::now();
		accept = 0;
		BOOST_FOREACH(Value const& v, lattice)
		{
			pgn = v.second;
			// cout << "Polygon " << k << endl;
			// Двигаем следующий полигон
			pgn_probe = MovePolygon(pgn, max_translation, max_angle, rnd, xmin, xmax, ymin, ymax);
			// Проверяем, что полигон pgn_probe ни с кем не пересекается (с учетом периодичности границ)
			lattice.remove(v);
			if (CheckPolygonIntersect(lattice, pgn_probe) == false)
			{
				Box_2 box;
				envelope(pgn_probe, box);
				lattice.insert(std::make_pair(box, pgn_probe));
				accept += 1;
			}
			else
				lattice.insert(v);
		};
		auto et = std::chrono::high_resolution_clock::now();
		cout << "Accept " << accept << " moves from " << lattice.size() << ". Step time - " << 
					std::chrono::duration_cast<std::chrono::milliseconds>(et-bt).count()/1000.0 << "second(s)." << endl;
		//SaveLattice(lattice, "step" + std::to_string(i));
	};
};

int main()
{
	// Радиус флека
	double rad_mul = 0.8;
	// Количество итераций при построении флека
	int iters_count = 2;
	// Вспомогательные переменные
	double sin30 = rad_mul*sin(M_PI/6);
	double sin60 = rad_mul*sin(M_PI/3);
	// Гексагон
	//Point hex_line[] = {Point(0, 0), Point(0, rad_mul), Point(-sin60, sin30), Point(-sin60, -sin30), Point(0, -rad_mul), Point(sin60, -sin30), 
						//Point(sin60, sin30), Point(0, rad_mul)};
	double hex_line[][2] = {{0, rad_mul}, {-sin60, sin30}, {-sin60, -sin30}, {0, -rad_mul}, {sin60, -sin30}, 
						{sin60, sin30}, {0, rad_mul}};


	// Построение флека
	std::vector<Point_2> base_polygon(0);
	for(int i = 0; i < 7; i++)
	{
		base_polygon.push_back(Point_2(hex_line[i][0], hex_line[i][1]));
	};
	for (int i = 0; i < iters_count; i++)
		base_polygon = HexFlakeStep(base_polygon);

	Polygon_2 hexaflake;
	Polygon_2 pgn1;
	Polygon_2 pgn2;
	append(hexaflake, base_polygon);
	// Размер решетки
	size_t nx = 128;
	size_t ny = 128;
	// Границы решетки
	double xmin = 0 - sin60;
	double xmax = (nx-1)*sqrt(3) + sqrt(3)/2 + sin60;
	double ymin = -2;
	double ymax = 1.5 + (ny-1)*3;
	// Решетка
	Lattice lattice;
	// Построение решетки
	Box_2 box1;
	Box_2 box2;
	for(int i=0; i < nx; i++) 
		for(int j=0; j < ny; j++)
	{
		strategy::transform::translate_transformer<double, 2, 2> t1(0 + i*sqrt(3), 0.5 + j*3);
		strategy::transform::translate_transformer<double, 2, 2> t2(sqrt(3)/2 + i*sqrt(3), -1 + j*3);
		transform(hexaflake, pgn1, t1);
		transform(hexaflake, pgn2, t2);
		envelope(pgn1, box1);
		envelope(pgn2, box2);
		lattice.insert(std::make_pair(box1, pgn1));
		lattice.insert(std::make_pair(box2, pgn2));
	}
	// Вывод флеков решетки
	//for(int i = 0; i < lattice.size(); i++)
		//cout << lattice[i] << endl;
	cout << "Run..." << endl;
	// Работа
	int steps_count = 100;
	Run(lattice, steps_count, nx, ny, xmin, xmax, ymin, ymax);
	cout << "Done!" << endl;
	// Сохранение флеков решетки
	//SaveLattice(lattice, "hexaflake");
  
    return 0;
}