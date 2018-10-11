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

#define pi 3.1415926535897932385

using namespace std;

ifstream in_file;
ofstream out_file;
int myid;// Номер потока
int numthreads;// Число потоков
int nx;// Общее число точек по пространственной оси
int nt;// Общее число точек по временной оси
int intnt;// Число точек для интегрирования по временной оси
double hx;// Шаг по пространственной оси
double ht;// Шаг по временной оси
double tmin;// Нижняя граница времени
double tmax;// Верхняя граница времени
double xmin;// Нижняя граница пространства
double xmax;// Верхняя граница пространства
int threads;// Количество потоков
double dmax;// Максимальная глобальная погрешность измерений на данном шаге
double dm;// Максимальная локальная погрешность измерений на данном шаге
double v;// Скорость первого импульса
double x0;// Координата центра первого импульса
double xb; // Координата дельта-барьера
double mu; // Мощность дельта-барьера
double l; // Обратная ширина дельта-барьера
double a;// Амплитуда ВЧ-поля
double w; // Частота ВЧ-поля
double** f;// Массив значений потенциала
int divx;// Делитель количества точек по оси координата при сохранении в файл
int divt;// Делитель количества точек по оси времени при сохранении в файл
double nu;// Коэффициент трения
string str;

// Начальный профиль
// Одиночный кинк
inline double f0_single_kink(const double x, const double t, const double v, const double x0)
{
	return 4*atan(exp((x-x0-v*t)/sqrt(1-v*v)));;
}

// Функция структурного возмущения
inline double delta_barrier(double x)
{
	double x2 = l*(x - xb)*l*(x - xb);
	return mu*l*exp(-x2)/sqrt(pi);
}

// Функция структурного возмущения
inline double right_part(double** f, int x, int t)
{
	double sum = 0.0;
	double A1, A;
	A = a*sin(w*ht*t);
	double t1b = t - intnt >= 0 ? t - intnt : 0;
	for(int t1 = 0; t1 <= t; t1++)
	{
		A1 = a*sin(w*ht*t1);
		sum += nu*ht*exp(-nu*(t-t1)*ht)*(sin(f[x][t1] + A1 - f[x][t] - A) + sin(f[x][t] + A));
	}
	return sum;
}

int main(int argc, char *argv[])
{
	//Параметры сетки и решения
	in_file.open(argv[1]);
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
	intnt = atoi(str.data());
	getline(in_file, str);
	getline(in_file, str);
	divx = atoi(str.data());
	getline(in_file, str);
	getline(in_file, str);
	divt = atoi(str.data());
	//Параметры импульсов
	getline(in_file, str);
	getline(in_file, str);
	v = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	x0 = atof(str.data());
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
	// Частота ВЧ поля
	getline(in_file, str);
	getline(in_file, str);
	w = atof(str.data());
	// Коэффициент трения
	getline(in_file, str);
	getline(in_file, str);
	nu = atof(str.data());
	in_file.close();
	l = 10;
	myid = 0;
	double begin = omp_get_wtime();
	cout << "Begin solve wave equation" << endl;
	//Начальные условия 
	hx = (xmax-xmin)/(nx-1);
	ht = hx/2;
	nt = (int)floor((tmax-tmin)/ht)+1; // Общее число узлов по времени
	// Потенциал решения
	f = (double**)calloc(nx, sizeof(double));//Число элементов по первому индексу
	for(int i = 0; i < nx; i++)
	{
		f[i] = (double*)calloc(nt, sizeof(double));//Число элементов по второму индексу
	};
	// Заполняем начальные условия
	for (int x = 0; x < nx; x++)
	{
		f[x][0] = f0_single_kink(xmin+x*hx, tmin, v, x0);
		f[x][1] = f0_single_kink(xmin+x*hx, tmin+ht, v, x0);
	}
	// Записываем начальный профиль в файл
	out_file.open(argv[2]);
	for(int x = 0; x < nx; x+=divx)
	{
		out_file << (f[x][1] - f[x][0])/ht << " ";
	};
	out_file << endl;
	out_file.close();
	omp_set_num_threads(threads);
	// Включаем время
	#pragma omp parallel private (myid) shared (f, threads, nx, divt, divx, nt, mu, xb, l, a, w)
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
		int k = 1;
		for(int t = 2; t < nt; t++)
		{
			for(int x = from_x; x < to_x; x++)
			{
				f[x][t] = (ht*ht)*(f[x-1][t-1]+f[x+1][t-1]-2*f[x][t-1])/(hx*hx)+2*f[x][t-1]-f[x][t-2]
						- ht*ht*((1 + delta_barrier(xmin + x*hx))*sin(f[x][t-1] + a*sin(w*ht*(t-1))) - right_part(f, x, t-1));
			};
			if(myid == 0)
				f[0][t] = (2*f[1][t]-f[2][t]);
			if(myid == threads-1)
				f[nx-1][t] = (2*f[nx-2][t]-f[nx-3][t]);
					if (myid == 0)
			if(myid == 0 && t % divt == 0)
			{
				cout << "Saving step #" << k << " from " << nt/divt  << endl;
				out_file.open(argv[2], std::ios_base::app);
				for(int x = 0; x < nx; x+=divx)
				{
					out_file << (f[x][t] - f[x][t - 1])/ht << " ";
				};
				out_file << endl;
				out_file.close();
				k += 1;
			};
			#pragma omp barrier
		};
	}
	double end = omp_get_wtime();
	cout << "Done solve wave equation" << endl;
	cout << "Computation took " << end - begin << " second(s)" << endl;
	//cin.get();
	return 0;
}

