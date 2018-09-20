#include <ctime>
#include <random>
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
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/algorithm/string.hpp>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
 
#define MAX_SOURCE_SIZE (0x1000000)

using namespace boost::geometry;
using namespace boost;
using std::cout; 
using std::endl;
BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian);

typedef model::d2::point_xy<float> Point_2;
typedef model::polygon<Point_2> Polygon_2;


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

void SaveLattice(std::string init_str, float * lattice, int count_poly, int coords_per_poly, std::string file_name)
{
	std::ofstream out_file;
	out_file.open (file_name + ".txt");
	out_file << init_str << endl;
	std::string str = "";
	for(int i = 0; i < count_poly; i++)
	{
        for(int k = 0; k < coords_per_poly; k+=2)
            out_file << std::to_string(lattice[i*coords_per_poly + k]) << " " << std::to_string(lattice[i*coords_per_poly + k + 1]) << " ";
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
inline void OutFlake(Polygon_2 pgn)
{
	std::string str = "";
	boost::geometry::for_each_point(pgn, [&str](Point_2 const& p) { str += std::to_string(get<0>(p)) + " " + std::to_string(get<1>(p)) + " "; });
	cout << str << endl;
}

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
	// Размер решетки
	int	nx =  atoi(argv[4]);
	int	ny = atoi(argv[4]);
	// Количество шагов по времени
	int steps_count = atoi(argv[5]);
	// Вспомогательные переменные
	float sin30 = rad_mul*sin(M_PI/6);
	float sin60 = rad_mul*sin(M_PI/3);
	float xmin, xmax, ymin, ymax;
	std::string init_str = "";
    int count_vertex_poly = 0;
	// Решетка
	std::vector<Polygon_2 > lattice(0);
	if (true)//(!file_exists("hexaflake_" + std::to_string(rad_mul) + "_" + std::to_string(iters_count)))
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
        count_vertex_poly = base_polygon.size();
		Polygon_2 hexaflake;
		Polygon_2 pgn1;
		Polygon_2 pgn2;
		append(hexaflake, base_polygon);
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
    // Пдготовка к работе (установка всех массивов для пересылки в вычислительное ядро)
    unsigned int count_poly = lattice.size();
    unsigned int coords_per_poly = count_vertex_poly * 2;
    unsigned int size = count_poly * coords_per_poly;
    float * lattice1 = (float*)malloc(sizeof(float) * size);
	float * lattice2 = (float*)malloc(sizeof(float) * size);
	float * centers = (float*)malloc(sizeof(float) * count_poly * 2);
	float * params = (float*)malloc(sizeof(float) * 6);
	unsigned int * seed = (unsigned int*)malloc(sizeof(unsigned int) * count_poly);
	unsigned int * accept = (unsigned int*)malloc(sizeof(unsigned int) * count_poly);
	params[0] = size;
	params[1] = count_poly;
	params[2] = xmin;
	params[3] = xmax;
	params[4] = ymin;
	params[5] = ymax;
	std::mt19937 gen(time(0));
	std::uniform_int_distribution<> uid(INT_MAX/2, INT_MAX);
	// Вспомогательные переменные
	Point_2 pc;
	int k = 0;
    for(int i = 0; i < count_poly; i++)
    {
		accept[i] = 0;
		seed[i] = uid(gen);
		centroid(lattice[i],pc);
		centers[2*i] = get<0>(pc);
		centers[2*i + 1] = get<1>(pc);
		k = 0;
        boost::geometry::for_each_point(lattice[i], 
            [&lattice1, &lattice2, &k, i, coords_per_poly](Point_2 const& p) { lattice1[i*coords_per_poly + k] = get<0>(p);
                                                                    lattice1[i*coords_per_poly + k + 1] = get<1>(p);
																	lattice2[i*coords_per_poly + k] = get<0>(p);
                                                                    lattice2[i*coords_per_poly + k + 1] = get<1>(p); k+=2; });
    }
	cout << "Run..." << endl;
	// Работа
	cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue command_queue = NULL;
    cl_mem lattice_dev = NULL;
	cl_mem probe_lattice_dev = NULL;
	cl_mem centers_dev = NULL;
	cl_mem params_dev = NULL;
	cl_mem seed_dev = NULL;
	cl_mem accept_dev = NULL;
    cl_program program = NULL;
    cl_kernel kernel = NULL;
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
    
    FILE *fp;
    char fileName[] = "./HexFlakeOCLKernel.c";
    char *source_str;
    size_t source_size;
    
    /* Load the source code containing the kernel*/
    fp = fopen(fileName, "r");
    if (!fp)
    {
        std::cout <<  "Failed to load kernel.\n";
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);
    
    /* Get Platform and Device Info */
    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    std::cout << "Platforms count - " << ret_num_platforms << std::endl;
    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
    std::cout << "Devices count - " << ret_num_devices << std::endl;
    /* Create OpenCL context */
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
    /* Create Command Queue */
   	command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
    /* Create Memory Buffers */
	lattice_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, size * sizeof(float), NULL, &ret);
	probe_lattice_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, size * sizeof(float), NULL, &ret);
	centers_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, 2 * count_poly * sizeof(float), NULL, &ret);
	params_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,  6 * sizeof(float), NULL, &ret);
	accept_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, count_poly * sizeof(int), NULL, &ret);
	seed_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, count_poly * sizeof(int), NULL, &ret);
	/* Copy data to the memory buffer */
    ret = clEnqueueWriteBuffer(command_queue, lattice_dev, CL_TRUE, 0, size * sizeof(float), lattice1, 0, NULL, NULL);
	ret = clEnqueueWriteBuffer(command_queue, probe_lattice_dev, CL_TRUE, 0, size * sizeof(float), lattice2, 0, NULL, NULL);
	ret = clEnqueueWriteBuffer(command_queue, centers_dev, CL_TRUE, 0, 2 * count_poly * sizeof(float), centers, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, params_dev, CL_TRUE, 0, 6 * sizeof(float), params, 0, NULL, NULL);
	ret = clEnqueueWriteBuffer(command_queue, accept_dev, CL_TRUE, 0, count_poly * sizeof(int), accept, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, seed_dev, CL_TRUE, 0, count_poly * sizeof(int), seed, 0, NULL, NULL);
    /* Create Kernel Program from the source */
    program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
    /* Build Kernel Program */
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    /* Create OpenCL Kernel */
    kernel = clCreateKernel(program, "Step", &ret);
    /* Set OpenCL Kernel Parameters */
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&lattice_dev);
	ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&probe_lattice_dev);
	ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&centers_dev);
	ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&params_dev);
	ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&accept_dev);
	ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&seed_dev);
    /* Execute OpenCL Kernel */
	size_t global_work_size = 16;//count_poly;
	int c = 0;
	while(c++ < steps_count*count_poly/global_work_size)
    	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	std::cout << "Run program result - " << ret << std::endl;
    /* Copy results from the memory buffer */
    ret = clEnqueueReadBuffer(command_queue, lattice_dev, CL_TRUE, 0, size * sizeof(float), lattice1, 0, NULL, NULL);
	ret = clEnqueueReadBuffer(command_queue, probe_lattice_dev, CL_TRUE, 0, size * sizeof(float), lattice2, 0, NULL, NULL);
	ret = clEnqueueReadBuffer(command_queue, centers_dev, CL_TRUE, 0, 2 * count_poly * sizeof(float), centers, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, params_dev, CL_TRUE, 0, 6 * sizeof(float), params, 0, NULL, NULL);
	ret = clEnqueueReadBuffer(command_queue, accept_dev, CL_TRUE, 0, count_poly * sizeof(int), accept, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, seed_dev, CL_TRUE, 0, count_poly * sizeof(int), seed, 0, NULL, NULL);    
    /* Finalization */
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(lattice_dev);
	ret = clReleaseMemObject(probe_lattice_dev);
	ret = clReleaseMemObject(centers_dev);
	ret = clReleaseMemObject(params_dev);
	ret = clReleaseMemObject(accept_dev);
	ret = clReleaseMemObject(seed_dev);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);
    
    free(source_str);
	cout << "Done!" << endl;
	// Сохранение флеков решетки
	init_str = SettingsToString(nx, ny, xmin, xmax, ymin, ymax);
	SaveLattice(init_str, lattice1, count_poly, coords_per_poly, "hexaflake_" + std::to_string(rad_mul) + "_" + std::to_string(iters_count));
  
    return 0;
}
