#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <chrono>
#include <omp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost::geometry;
using namespace boost;
using std::cout; 
using std::endl;
BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian);

typedef model::d2::point_xy<float> Point_2;
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
	float Nextfloat()
	{
		return ((float)xor128()) / ULONG_MAX;
	};
};

bool file_exists (const std::string name) {
	return ( access( (name + ".txt").c_str(), F_OK ) != -1 );
}

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
        float dpx = (get<0>(p2) - get<0>(p1))/3.0;
        float dpy = (get<1>(p2) - get<1>(p1))/3.0;
		new_line.push_back(p1);
        p3 = Point_2(get<0>(p1) + dpx, get<1>(p1) + dpy);
        p4 = Point_2(get<0>(p1) + 2*dpx, get<1>(p1) + 2*dpy);
        float alpha = M_PI/3.0;
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


// Функция сохранения флеков решетки lattice в файл file_name
void SaveLattice(std::string init_str, std::vector<Polygon_2> &lattice, std::string file_name)
{
	std::ofstream out_file;
	out_file.open (file_name + ".txt");
	out_file << init_str << endl;
	std::string str = "";
	for(int i = 0; i < lattice.size(); i++)
	{
		boost::geometry::for_each_point(lattice[i], [&out_file](Point_2 const& p) { out_file << std::to_string(get<0>(p)) + " " + std::to_string(get<1>(p)) << " "; });
		out_file << endl;
	}
	out_file.close();
}

// Функция загрузки флеков решетки lattice из файла file_name
std::vector<Polygon_2> LoadLattice(std::string file_name, std::string &init_str)
{
	std::vector<Polygon_2> lattice(0);
	std::ifstream in_file;
	in_file.open (file_name + ".txt");
	std::string str = "";
	getline(in_file, init_str);
	cout << init_str << endl;
	while(getline(in_file, str)){
		trim(str);
		std::vector<std::string> results(0);
		boost::algorithm::split(results, str, is_any_of(" "));
		std::vector<Point_2> base_polygon(0);
		for(int i = 0; i < results.size(); i += 2)
		{
			base_polygon.push_back(Point_2(atof(results[i].c_str()), atof(results[i+1].c_str())));
		};
		Polygon_2 hexaflake;
		append(hexaflake, base_polygon);
		lattice.push_back(hexaflake);
	}
	in_file.close();
	return lattice;
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
						int pgn_index, float max_translation, float max_angle, RandomGenerator &rnd, 
						float xmin, float xmax, float ymin, float ymax)
{
	Polygon_2 ret, ret1;
	if( rnd.Nextfloat() > 0.5)
	{
		float dx = -get<0>(centers[pgn_index]);
		float dy = -get<1>(centers[pgn_index]);
		float angle = max_angle*(1 - 2*rnd.Nextfloat());
	    strategy::transform::translate_transformer<float, 2, 2> t1(dx, dy);
		transform(lattice[pgn_index], ret, t1);
	    strategy::transform::rotate_transformer<boost::geometry::radian, float, 2, 2> r(angle);
		transform(ret, ret1, r);
	    strategy::transform::translate_transformer<float, 2, 2> t2(-dx, -dy);
		transform(ret1, ret, t2);
	}
	else{
		float x = get<0>(centers[pgn_index]);
		float y = get<1>(centers[pgn_index]);
		float dx = max_translation*(1 - 2*rnd.Nextfloat());
		float dy = max_translation*(1 - 2*rnd.Nextfloat());
		if (x + dx > xmax)
			dx = dx - (xmax - xmin);
		else if (x + dx < xmin)
			dx = dx + (xmax - xmin);
		if (y + dy > ymax)
			dy = dy - (ymax - ymin);
		else if (y + dy < ymin)
			dy = dy + (ymax - ymin);
	    strategy::transform::translate_transformer<float, 2, 2> t(dx, dy);
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
	float dist = 0;
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
		int steps_count, int nx, int ny, float xmin, float xmax, float ymin, float ymax, int count_threads)
{
	RandomGenerator rnd;
	Polygon_2 pgn_probe;
	std::vector<Point_2 > centers(lattice.size());
	for(int i=0; i < lattice.size(); i++)
	{
		centroid(lattice[i], centers[i]);
	};
	// Максимальные смещения
	float max_translation = 0.1;
	float max_angle = 0.1;
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

std::string SettingsToString(int nx, int ny, float xmin, float xmax, float ymin, float ymax)
{
	return std::to_string(nx) + " " +std::to_string(ny) + " " +std::to_string(xmin) + " " + 
		std::to_string(xmax) + " " + std::to_string(ymin) + " " + std::to_string(ymax); 
}

void SettingsFromString(std::string str, int& nx, int& ny, float& xmin, float& xmax, float& ymin, float& ymax)
{
	std::vector<std::string> results;
	boost::algorithm::split(results, str, is_any_of(" "));
	nx = atoi(results[0].c_str());
	ny = atoi(results[1].c_str());
	xmin = atof(results[2].c_str());
	xmax = atof(results[3].c_str());
	ymin = atof(results[4].c_str());
	ymax = atof(results[5].c_str());
}

int main(int argc, char* argv[])
{
	// Радиус флека
	float rad_mul = atof(argv[1]);
	// Количество итераций при построении флекса
	int iters_count = atoi(argv[2]);
	// Устанавливаем количество потоков
	int count_threads = atoi(argv[3]);
	// Вспомогательные переменные
	float sin30 = rad_mul*sin(M_PI/6);
	float sin60 = rad_mul*sin(M_PI/3);
	int nx, ny;
	float xmin, xmax, ymin, ymax;
	std::string init_str = "";
	// Решетка
	std::vector<Polygon_2 > lattice(0);
	if (!file_exists("hexaflake_" + std::to_string(rad_mul) + "_" + std::to_string(iters_count)))
	{
		std::cout << "First run. Construct hexflake..." << endl;
		// Гексагон
		//Point hex_line[] = {Point(0, 0), Point(0, rad_mul), Point(-sin60, sin30), Point(-sin60, -sin30), Point(0, -rad_mul), Point(sin60, -sin30), 
							//Point(sin60, sin30), Point(0, rad_mul)};
		float hex_line[][2] = {{0, rad_mul}, {-sin60, sin30}, {-sin60, -sin30}, {0, -rad_mul}, {sin60, -sin30}, 
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
		nx =  atoi(argv[4]);
		ny = atoi(argv[4]);
		// Границы решетки
		xmin = 0 - sin60;
		xmax = (nx-1)*sqrt(3) + sqrt(3)/2 + sin60;
		ymin = -2;
		ymax = 1.5 + (ny-1)*3;
		// Построение решетки
		for(int i=0; i < nx; i++) 
			for(int j=0; j < ny; j++)
		{
			strategy::transform::translate_transformer<float, 2, 2> t1(0 + i*sqrt(3), 0.5 + j*3);
			strategy::transform::translate_transformer<float, 2, 2> t2(sqrt(3)/2 + i*sqrt(3), -1 + j*3);
			transform(hexaflake, pgn1, t1);
			transform(hexaflake, pgn2, t2);
			lattice.push_back(pgn1);
			lattice.push_back(pgn2);
		}
	}
	else
	{
		std::cout << "Load hexflake from file..." << endl;
		lattice = LoadLattice("hexaflake_" + std::to_string(rad_mul) + "_" + std::to_string(iters_count), init_str);
		SettingsFromString(init_str, nx, ny, xmin, xmax, ymin, ymax);
	}
	cout << "Run..." << endl;
	// Работа
	int steps_count = atoi(argv[5]);
	Run(lattice, steps_count, nx, ny, xmin, xmax, ymin, ymax, count_threads);
	cout << "Done!" << endl;
	// Сохранение флеков решетки
	init_str = SettingsToString(nx, ny, xmin, xmax, ymin, ymax);
	SaveLattice(init_str, lattice, "hexaflake_" + std::to_string(rad_mul) + "_" + std::to_string(iters_count));
  
    return 0;
}
