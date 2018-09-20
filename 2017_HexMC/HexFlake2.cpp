#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/squared_distance_2.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Aff_transformation_2<K> Transformation; 
typedef CGAL::Vector_2<K> Vector; 
typedef CGAL::Direction_2<K> Direction; 


using std::cout; 
using std::endl;

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
Polygon_2 HexFlakeStep(Polygon_2 &old_line)
{
    std::vector<Point> new_line(0);
	size_t old_line_len = old_line.size();
	Point p, p1, p2, p3, p4, p5;
	new_line.push_back(old_line[0]);
    for( size_t i = 1; i < old_line_len-1; i++){
        p1 = old_line[i];
        p2 = old_line[i+1];
        double dpx = (p2[0] - p1[0])/3.0;
        double dpy = (p2[1] - p1[1])/3.0;
        new_line.push_back(p1);
        p3 = Point(p1[0] + dpx, p1[1] + dpy);
        p4 = Point(p1[0] + 2*dpx, p1[1] + 2*dpy);
        double alpha = M_PI/3.0;
        p = Point((p2[0] - p1[0])/3.0, (p2[1] - p1[1])/3.0);
        p5 = Point(p[0]*cos(alpha) - p[1]*sin(alpha), p[1]*cos(alpha) + p[0]*sin(alpha));
        p5 = Point(p5[0] + p3[0], p5[1] + p3[1]);
        new_line.push_back(p3);
        new_line.push_back(p5);
        new_line.push_back(p4);
	};
    new_line.push_back(old_line[old_line_len-1]);
	Point *mas = new Point[new_line.size()];
	for(size_t i = 0; i < new_line.size(); i++)
		mas[i] = new_line[i];
	Polygon_2 ret(mas, mas + new_line.size());
    return  ret;
}

// Функция сохранения флеков решетки lattice в файл "hexaflake.txt"
void SaveLattice(std::vector<Polygon_2> lattice, std::string file_name)
{
	std::ofstream out_file;
	out_file.open (file_name + ".txt");
	for(int i = 0; i < lattice.size(); i++)
		out_file << lattice[i] << endl;
	out_file.close();
}

// Двигаем полигон с индексом pgn_index
Polygon_2 MovePolygon(std::vector<Polygon_2> &lattice, int pgn_index, double max_translation, double max_angle, RandomGenerator &rnd, 
						double xmin, double xmax, double ymin, double ymax)
{
	Polygon_2 ret;
	if( rnd.NextDouble() > 0.5)
	{
		double angle = max_angle*(1 - 2*rnd.NextDouble());
	    Transformation r(CGAL::ROTATION, sin(angle), cos(angle));
		ret = CGAL::transform(r, lattice[pgn_index]);
	}
	else{
		double x = lattice[pgn_index].vertex(0).x();
		double y = lattice[pgn_index].vertex(0).y();
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
	    Transformation t(CGAL::TRANSLATION, Vector(dx, dy));
		ret = CGAL::transform(t, lattice[pgn_index]);		
	}
	return ret;
}

// Проверяем, что полигон pgn_probe, полученный из полигона с индексом pgn_index, никого не перекрывает
bool CheckPolygonIntersect(std::vector<Polygon_2> &lattice, Polygon_2 &pgn_probe, int pgn_index)
{
	bool flag = false;
	double dist = 0;
	for(int i=0; i < lattice.size(); i++)
	{
		dist = CGAL::squared_distance(lattice[i].vertex(0), pgn_probe.vertex(0));
		//cout << dist << endl;
		if ( i == pgn_index || dist > 4 )
		{
			continue;
		}
		else
		{
			//cout << pgn_probe << endl;
			//cout << lattice[i] << endl;
			
			flag = CGAL::do_intersect(lattice[i], pgn_probe);
			//cout << flag << endl;
			if ( flag == true )
				break;
		};
	};			
	return flag;
}

void Run(std::vector<Polygon_2> &lattice, int steps_count, int nx, int ny, double xmin, double xmax, double ymin, double ymax)
{
	RandomGenerator rnd;
	Polygon_2 pgn_probe;
	// Максимальные смещения
	double max_translation = 0.1;
	double max_angle = 0.1;
	for(int i = 0; i < steps_count; i++)
	{
		cout << "Time step " << i << endl;
		for(int k = 0; k < lattice.size(); k++)
		{
			// cout << "Polygon " << k << endl;
			// Двигаем полигон с индексом k
			pgn_probe = MovePolygon(lattice, k, max_translation, max_angle, rnd, xmin, xmax, ymin, ymax);
			// Проверяем, что полигон pgn_probe ни с кем не пересекается (с учетом периодичности границ)
			if (CheckPolygonIntersect(lattice, pgn_probe, k) == false)
			{
				lattice[k] = pgn_probe;
			};
		};
		SaveLattice(lattice, std::to_string(i));
	};
};

int main()
{
	// Радиус флекса
	double rad_mul = 0.8;
	// Количество итераций при построении флекса
	int iters_count = 0;
	// Вспомогательные переменные
	double sin30 = rad_mul*sin(M_PI/6);
	double sin60 = rad_mul*sin(M_PI/3);
	// Гексагон
	//Point hex_line[] = {Point(0, 0), Point(0, rad_mul), Point(-sin60, sin30), Point(-sin60, -sin30), Point(0, -rad_mul), Point(sin60, -sin30), 
						//Point(sin60, sin30), Point(0, rad_mul)};
	Point hex_line[] = {Point(0, rad_mul), Point(-sin60, sin30), Point(-sin60, -sin30), Point(0, -rad_mul), Point(sin60, -sin30), 
						Point(sin60, sin30)};

	Polygon_2 hexagone(hex_line, hex_line + 6);
	// Построение флека
	Polygon_2 base_polygon = hexagone;
	for (int i = 0; i < iters_count; i++)
		base_polygon = HexFlakeStep(base_polygon);
	// Размер решетки
	size_t nx = 10;
	size_t ny = 10;
	// Границы решетки
	double xmin = 0 - sin60;
	double xmax = (nx-1)*sqrt(3) + sqrt(3)/2 + sin60;
	double ymin = -2;
	double ymax = 1.5 + (ny-1)*3;
	// Решетка
	std::vector<Polygon_2> lattice(0);
	// Построение решетки
	for(int i=0; i < nx; i++) 
		for(int j=0; j < ny; j++)
	{
		Transformation t1(CGAL::TRANSLATION, Vector(0 + i*sqrt(3), 0.5 + j*3));
		Transformation t2(CGAL::TRANSLATION, Vector(sqrt(3)/2 + i*sqrt(3), -1 + j*3));
		Polygon_2 pgn1 = CGAL::transform( t1, base_polygon);
		Polygon_2 pgn2 = CGAL::transform( t2, base_polygon);
		lattice.push_back(pgn1);
		lattice.push_back(pgn2);
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
	SaveLattice(lattice, "hexaflake");
  
    return 0;
}