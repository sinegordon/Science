#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <omp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>

using namespace boost::geometry;
using std::cout; 
using std::endl;
BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian);

typedef model::d2::point_xy<double> Point_2;
typedef model::polygon<Point_2> Polygon_2;


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
void SaveLattice(std::vector<Polygon_2> &lattice, std::string file_name)
{

	std::ofstream out_file;
	out_file.open (file_name + ".txt");
	std::string str = "";
	for(int i = 0; i < lattice.size(); i++)
	{
		boost::geometry::for_each_point(lattice[i], [&out_file](Point_2 const& p) { out_file << std::to_string(get<0>(p)) + " " + std::to_string(get<1>(p)) << " "; });
		out_file << endl;
	}
	out_file.close();
}

// Функция вывода флека  в стандартный поток
void OutFlake(Polygon_2 pgn)
{
	std::string str = "";
	boost::geometry::for_each_point(pgn, [&str](Point_2 const& p) { str += std::to_string(get<0>(p)) + " " + std::to_string(get<1>(p)) + " "; });
	cout << str << endl;
}

// Двигаем полигон с индексом pgn_index
Polygon_2 MovePolygon(std::vector<Polygon_2> &lattice, 
						std::vector<Point_2> &centers,
						int pgn_index, double max_translation, double max_angle, RandomGenerator &rnd, 
						double xmin, double xmax, double ymin, double ymax)
{
	Polygon_2 ret, ret1;
	if( rnd.NextDouble() > 0.5)
	{
		double dx = -get<0>(centers[pgn_index]);
		double dy = -get<1>(centers[pgn_index]);
		double angle = max_angle*(1 - 2*rnd.NextDouble());
	    strategy::transform::translate_transformer<double, 2, 2> t1(dx, dy);
		transform(lattice[pgn_index], ret, t1);
	    strategy::transform::rotate_transformer<boost::geometry::radian, double, 2, 2> r(angle);
		transform(ret, ret1, r);
	    strategy::transform::translate_transformer<double, 2, 2> t2(-dx, -dy);
		transform(ret1, ret, t2);
	}
	else{
		double x = get<0>(centers[pgn_index]);
		double y = get<1>(centers[pgn_index]);
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
		transform(lattice[pgn_index], ret, t);
	}
	return ret;
}

// Проверяем, что полигон pgn_probe, полученный из полигона с индексом pgn_index, никого не перекрывает
bool CheckPolygonIntersect(std::vector<Polygon_2> &lattice,
						   std::vector<Point_2> &centers,
	 					   Polygon_2 &pgn_probe, int pgn_index)
{
	bool flag = false;
	double dist = 0;
	Point_2 p;
	centroid(pgn_probe, p);
	for(int i=0; i < lattice.size(); i++)
	{
		dist = comparable_distance(p, centers[i]);
		//cout << dist << endl;
		if ( i == pgn_index || dist > 4 )
		{
			continue;
		}
		else
		{
			//cout << pgn_probe << endl;
			//cout << lattice[i] << endl;
			
			flag = overlaps(lattice[i], pgn_probe);
			//cout << flag << endl;
			if ( flag == true )
				break;
		};
	};			
	return flag;
}

void Run(std::vector<Polygon_2> &lattice, 
		int steps_count, int nx, int ny, double xmin, double xmax, double ymin, double ymax, int count_threads)
{
	RandomGenerator rnd;
	Polygon_2 pgn_probe;
	std::vector<Point_2 > centers(lattice.size());
	for(int i=0; i < lattice.size(); i++)
	{
		centroid(lattice[i], centers[i]);
	};
	// Максимальные смещения
	double max_translation = 0.1;
	double max_angle = 0.1;
	int accept_global = 0, accept_local = 0;
	omp_set_num_threads(count_threads);
	std::vector<Polygon_2> pgn_probes(count_threads);
	for(int i = 0; i < steps_count; i++)
	{
		cout << "Time step " << i << endl;
		auto bt = std::chrono::high_resolution_clock::now();
		accept_global = 0;
		//bool flag = false;
		#pragma omp parallel private(pgn_probe, accept_local) shared(pgn_probes, accept_global)
		{
			int size = lattice.size();
			int num = omp_get_thread_num();
			int low = num*size/count_threads;
			int high = (num + 1)*size/count_threads;
			accept_local = 0;
			for (int k = low; k < high; k++)
			{
				// cout << "Polygon " << k << endl;
				// Двигаем полигон с индексом k пока он ни с кем не пересекается
				pgn_probes[num] = MovePolygon(lattice, centers, k, max_translation, max_angle, rnd, xmin, xmax, ymin, ymax);
				/*
				do{
					pgn_probes[num] = MovePolygon(lattice, centers, k, max_translation, max_angle, rnd, xmin, xmax, ymin, ymax);
					#pragma omp barrier
					if(num == 0)
					{
						flag = false;					
						for(int t = 0; t < count_threads; t++)
							for(int l=0; l < count_threads; l++)
								flag = flag || overlaps(pgn_probes[t], pgn_probes[l]);
					
					};
					#pragma omp barrier
				} while (flag == true);
				*/
				#pragma omp barrier
				// Проверяем, что полигон pgn_probe ни с кем не пересекается (с учетом периодичности границ)
				if (CheckPolygonIntersect(lattice, centers, pgn_probes[num], k) == false)
				{
					lattice[k] = pgn_probes[num];
					boost::geometry::centroid(lattice[k], centers[k]);
					accept_local += 1;
				};
			};
			#pragma omp atomic
			accept_global += accept_local;
		};
		auto et = std::chrono::high_resolution_clock::now();
		cout << "Accept " << accept_global << " moves from " << lattice.size() << ". Step time - " << 
					std::chrono::duration_cast<std::chrono::milliseconds>(et-bt).count()/1000.0 << "second(s)." << endl;
		//SaveLattice(lattice, "step" + std::to_string(i));
	};
};

int main(int argc, char* argv[])
{
	// Радиус флека
	double rad_mul = atof(argv[1]);
	// Количество итераций при построении флекса
	int iters_count = atoi(argv[2]);
	// Устанавливаем количество потоков
	int count_threads = atoi(argv[3]);
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
	size_t nx =  atoi(argv[4]);
	size_t ny = atoi(argv[4]);
	// Границы решетки
	double xmin = 0 - sin60;
	double xmax = (nx-1)*sqrt(3) + sqrt(3)/2 + sin60;
	double ymin = -2;
	double ymax = 1.5 + (ny-1)*3;
	// Решетка
	std::vector<Polygon_2 > lattice(0);
	// Построение решетки
	for(int i=0; i < nx; i++) 
		for(int j=0; j < ny; j++)
	{
		strategy::transform::translate_transformer<double, 2, 2> t1(0 + i*sqrt(3), 0.5 + j*3);
		strategy::transform::translate_transformer<double, 2, 2> t2(sqrt(3)/2 + i*sqrt(3), -1 + j*3);
		transform(hexaflake, pgn1, t1);
		transform(hexaflake, pgn2, t2);
		lattice.push_back(pgn1);
		lattice.push_back(pgn2);
	}
	// Вывод флеков решетки
	//for(int i = 0; i < lattice.size(); i++)
		//cout << lattice[i] << endl;
	cout << "Run..." << endl;
	// Работа
	int steps_count = atoi(argv[5]);
	Run(lattice, steps_count, nx, ny, xmin, xmax, ymin, ymax, count_threads);
	cout << "Done!" << endl;
	// Сохранение флеков решетки
	SaveLattice(lattice, "hexaflake_" + std::to_string(rad_mul) + "_" + std::to_string(iters_count));
  
    return 0;
}