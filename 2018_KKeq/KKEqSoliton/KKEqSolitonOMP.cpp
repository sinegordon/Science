// Soliton.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <limits>
#include "alglib/src/interpolation.h"

#define pi 3.1415926535897932385

using namespace std;
using namespace alglib;

ifstream in_file;
ofstream out_file;
int myid;// Номер потока
int numthreads;// Число потоков
int nx;// Общее число точек по пространственной оси
int nt;// Общее число точек по временной оси
int masnt;// Число точек в массиве по временной оси
double hx;// Шаг по пространственной оси
double ht;// Шаг по временной оси
double tmin;// Нижняя граница времени
double tmax;// Верхняя граница времени
double xmin;// Нижняя граница пространства
double xmax;// Верхняя граница пространства
int threads;// Количество потоков
double dmax;// Максимальная глобальная погрешность измерений на данном шаге
double dm;// Максимальная локальная погрешность измерений на данном шаге
double v1;// Скорость первого импульса
double v2;// Скорость второго импульса
double x10;// Координата центра первого импульса
double x20;// Координата центра второго импульса
double xb; // Координата дельта-барьера
double mu; // Мощность дельта-барьера
double l; // Обратная ширина дельта-барьера
double a;// Амплитуда ВЧ-поля
double b;// Параметр b (отношение энергий)
double** f;// Массив значений потенциала
double* xmas;// Массив точек по оси абсцисс в файле с таблично заданным кинком
double* ymas;// Массив точек по оси ординат в файле с таблично заданным кинком
double xmin_kink;// Минимальное значение интервала интерполяции
double xmax_kink;// Максимальное значение интервала интерполяции
double* xmas1;// Массив точек по оси абсцисс в файле с таблично заданным током
double* ymas1;// Массив точек по оси ординат в файле с таблично заданным током
double fmin_kink;// Минимальное значение потенциала при интерполяции
double fmax_kink;// Максимальное значение потенциала при интерполяции
int size_kink;// Количество точек в таблице задания кинка
int divx;// Делитель количества точек по оси координата при сохранении в файл
int divt;// Делитель количества точек по оси времени при сохранении в файл
string str;
double * alpha_hf;// Вспомогательная переменная для потенциала при усреднении по периоду ВЧ поля
double F0; //  F(0,a,b) в подкоренном выражении интеграла для усредненного кинка
// Пременные для построения сплайна кинка
real_1d_array xa;
real_1d_array ya;
barycentricinterpolant p;
spline1dinterpolant s;
// Пременные для построения сплайна тока
real_1d_array xa1;
real_1d_array ya1;
barycentricinterpolant p1;
spline1dinterpolant s1;

// Подынтегральная функция при усредении по периоду ВЧ поля
void int_function_HF(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = sqrt(1 + b*b*(1 - cos(alpha_hf[myid] + a*sin(x))));
}

// Интеграл, определяющий значение функции в подкоренном выражении интеграла для кинка (F(alpha, a, b) в статье)
double F(double alpha)
{
    autogkstate s;
    double v = 0.0;
    autogkreport rep;
    autogksmooth(0, 2*pi, s);
	alpha_hf[myid] = alpha;
    alglib::autogkintegrate(s, int_function_HF);
    autogkresults(s, v, rep);
	return v/2.0/pi;
}

// Подынтегральная функция при неявном определении кинка усредненного по периоду ВЧ поля
void int_function_av_kink(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(F(x) - F0);
}

// Подынтегральная функция при неявном определении кинка
void int_function_kink(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(sqrt(1 + b*b*(1 - cos(x))) - 1);
}

// Интеграл, определяющий значение кинка в неявном виде
double Fun(double alpha)
{
    autogkstate s;
    double v = 0.0;
    autogkreport rep;
    autogksmooth(pi, alpha, s);
    alglib::autogkintegrate(s, int_function_kink);
    autogkresults(s, v, rep);
	return v;
}

// Интеграл, определяющий значение усредененного кинка в неявном виде
double FunHF(double alpha)
{
    autogkstate s;
    double v = 0.0;
    autogkreport rep;
    autogksmooth(pi, alpha, s);
    alglib::autogkintegrate(s, int_function_av_kink);
    autogkresults(s, v, rep);
	return v;
}

// Начальный профиль
// Одиночный кинк
double f0_single_kink(const double x, const double t, const double v1, const double x10, const double b)
{

	double f1 = 0.0;
	double arg1 = (x-x10-v1*t)/sqrt(1-v1*v1);

	if(arg1 <= xmin_kink)
	{
		f1 = spline1dcalc(s, xmin_kink);
	}
	else
	{
		if (arg1 >= xmax_kink)
		{
			f1 = spline1dcalc(s, xmax_kink);
		}
		else
		{
			f1 = spline1dcalc(s, arg1);
		}
	};
	return f1;
}

// Двойной кинк
double f0_double_kink(const double x, const double t, const double v1, const double v2, const double x10, const double x20, const double b)
{

	double f1 = 0.0, f2 = 0.0;
	double arg2 = (x-x20-v2*t)/sqrt(1-v2*v2);
	double arg1 = (x-x10-v1*t)/sqrt(1-v1*v1);

	if(arg1 <= xmin_kink)
	{
		f1 = 0.0;
	}
	else
	{
		if (arg1 >= xmax_kink)
		{
			f1 = 2*pi;
		}
		else
		{
			f1 = spline1dcalc(s, arg1);
		}
	};
	if(arg2 <= xmin_kink)
	{
		f2 = 0.0;
	}
	else
	{
		if (arg2 >= xmax_kink)
		{
			f2 = 2*pi;
		}
		else
		{
			f2 = spline1dcalc(s, arg2);
		}
	};
	return f1 + f2;
}

// Функция структурного возмущения
inline double delta_barrier(double x)
{
	return mu*l*exp(-sqr(l*(x - xb)))/sqrt(pi);
}


// Подынтегральная функция для тока усредненного по периоду ВЧ поля
void int_function_av_current(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
	y = b*b*sin(alpha_hf[myid] + a*sin(x))/sqrt(1 + b*b*(1 - cos(alpha_hf[myid] + a*sin(x))));
};

// Функция тока с усреднением по ВЧ полю
double currentHF(double f)
{
	autogkstate s;
	double v = 0.0;
	autogkreport rep;
	autogksmooth(0, 2*pi, s);
	alpha_hf[myid] = f;
	alglib::autogkintegrate(s, int_function_av_current);
	autogkresults(s, v, rep);
	return v/2.0/pi;
};

double current(double f)
{
	return b*b*sin(f)/sqrt(1 + b*b*(1 - cos(f)));
};

int main(int argc, char *argv[])
{
	//Параметры сетки и решения
	in_file.open("inOMP.txt");
	getline(in_file, str);
	getline(in_file, str);
	tmin = atof(str.data());	
	getline(in_file, str);
	tmax = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	xmin = atof(str.data());	
	getline(in_file, str);
	xmax = atof(str.data());
	getline(in_file, str);
	nx = atof(str.data());
	getline(in_file, str);
	getline(in_file, str);
	threads = atoi(str.data());
	getline(in_file, str);
	getline(in_file, str);
	b = atof(str.data());
	getline(in_file, str);
	getline(in_file, str);
	masnt = atoi(str.data());
	getline(in_file, str);
	getline(in_file, str);
	divx = atoi(str.data());
	getline(in_file, str);
	getline(in_file, str);
	divt = atoi(str.data());
	//Параметры импульсов
	getline(in_file, str);
	getline(in_file, str);
	v1 = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	x10 = atof(str.data());
	getline(in_file, str);
	getline(in_file, str);
	v2 = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	x20 = atof(str.data());
	// Параметры дельта-барьера
	getline(in_file, str);
	getline(in_file, str);
	xb = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	mu = atof(str.data());
	// Амплитуда ВЧ поля
	getline(in_file, str);
	getline(in_file, str);
	a = atof(str.data());	
	in_file.close();
	l = 10;
	cout << "Begin interpolation kink and current" << endl;
	alpha_hf = new double[threads];
	myid = 0;
	F0 = F(0);
	double begin = omp_get_wtime();
	// Построение таблично заданного кинка
	// Строим неравномерную таблицу с помощью неявно заданной формы кинка
	vector<double> xmas_temp; 
	vector<double> ymas_temp;
	int k = 0;
	for(double alpha = 0.000001; alpha < 2*pi; alpha += 0.0005)
	{
		xmas_temp.push_back(FunHF(alpha)/2);
		ymas_temp.push_back(alpha);
		if(k % 100 == 0)
			cout << "Kink alpha = " << alpha << endl;
		k += 1;
	}
	xmax_kink = max(xmas_temp[0], xmas_temp[xmas_temp.size()-1]);
	xmin_kink = min(xmas_temp[0], xmas_temp[xmas_temp.size()-1]);
	size_kink = xmas_temp.size();
	xmas = new double[size_kink];
	ymas = new double[size_kink];
	xmas1 = new double[size_kink];
	ymas1 = new double[size_kink];
	k = 0;
	for(int i = 0; i < size_kink; i++)
	{
		xmas[i] = xmas_temp[i];
		ymas[i] = ymas_temp[i];
		xmas1[i] = ymas_temp[i];
		ymas1[i] = currentHF(xmas1[i]);
		if(k % 100 == 0)
			cout << "Current alpha = " << xmas1[i] << endl;
		k += 1;
	}
	// Строим сплайн, соответствующий таблично заданному кинку
	xa.setcontent(size_kink, xmas);
	ya.setcontent(size_kink, ymas);
	spline1dbuildlinear(xa, ya, s);
	// Строим сплайн, соответствующий таблично заданному току
	xa1.setcontent(size_kink, xmas1);
	ya1.setcontent(size_kink, ymas1);
	spline1dbuildlinear(xa1, ya1, s1);
	cout << "Done interpolation" << endl;
	cout << "Begin solve wave equation" << endl;
	//Начальные условия 
	hx = (xmax-xmin)/(nx-1);
	ht = hx/2;
	nt = (int)floor((tmax-tmin)/ht)+1; // Общее число узлов по времени
	f = (double**)calloc(nx, sizeof(double));//Число элементов по первому индексу
	for(int i = 0; i < nx; i++)
	{
		f[i] = (double*)calloc(masnt, sizeof(double));//Число элементов по второму индексу
	};
	// Заполняем начальные условия
	// Если одиночный кинк
	if (v2 == 0.0)
		for (int x = 0; x < nx; x++)
		{
			f[x][0] = f0_single_kink(xmin+x*hx, tmin, v1, x10, b);
			f[x][1] = f0_single_kink(xmin+x*hx, tmin+ht, v1, x10, b);
		}
	else
		// Если кинк двойной
		for (int x = 0; x < nx; x++)
		{
			f[x][0] = f0_double_kink(xmin+x*hx, tmin, v1, v2, x10, x20, b);
			f[x][1] = f0_double_kink(xmin+x*hx, tmin+ht, v1, v2, x10, x20, b);
		};
	omp_set_num_threads(threads);
	// Включаем время
	#pragma omp parallel private (myid) shared (f, threads, masnt, nx, divt, divx, nt, mu, xb, l, alpha_hf)
	{
		myid = omp_get_thread_num();
		int from_x, to_x;
		if(myid > 0)
			from_x = myid*nx/threads;
		else
			from_x = 1;
		if(myid < threads-1)
			to_x = (myid+1)*nx/threads;
		else
			to_x = nx-1;
		for(int k = 0; k < nt / masnt; k++)
		{
			for(int t = 2; t < masnt; t++)
			{
				for(int x = from_x; x < to_x; x++)
				{
					f[x][t] = (ht*ht)*(f[x-1][t-1]+f[x+1][t-1]-2*f[x][t-1])/(hx*hx)+2*f[x][t-1]-f[x][t-2]
							- (1 + delta_barrier(xmin + x*hx))*ht*ht*spline1dcalc(s1, f[x][t-1]);//spline1dcalc - ток
				};
				if(myid == 0)
					f[0][t] = (2*f[1][t]-f[2][t]);
				if(myid == threads-1)
					f[nx-1][t] = (2*f[nx-2][t]-f[nx-3][t]);
				#pragma omp barrier
			};
			if (myid == 0)
			{
				cout << "Saving iteration #" << k << " from " << nt / masnt << endl;
				if(k == 0)
					out_file.open("out_f.sgo");
				else
					out_file.open("out_f.sgo", std::ios_base::app);
				for(int t = 0; t < masnt; t+=divt)
				{
					for(int x = 0; x < nx - 1; x+=divx)
					{
						out_file << f[x][t] << " ";
					};
					out_file << f[nx - 1][t] << endl;
				};
				out_file.close();
				for (int x = 0; x < nx; x++)
				{
					f[x][0] = f[x][masnt-2];
					f[x][1] = f[x][masnt-1];
				};
			}
			#pragma omp barrier
		};
	}
	double end = omp_get_wtime();
	cout << "Done solve wave equation" << endl;
	cout << "Computation took " << end - begin << " second(s)" << endl;
	cin.get();
	return 0;
}

