// Soliton.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <limits>
#include "alglib\src\interpolation.h"

#define pi 3.1415926535897932385

using namespace std;
using namespace alglib;

ifstream in_file;
ofstream out_file;
int myid;// ����� ������
int numthreads;// ����� �������
int nx;// ����� ����� ����� �� ���������������� ���
int nt;// ����� ����� ����� �� ��������� ���
int masnt;// ����� ����� � ������� �� ��������� ���
double hx;// ��� �� ���������������� ���
double ht;// ��� �� ��������� ���
double tmin;// ������ ������� �������
double tmax;// ������� ������� �������
double xmin;// ������ ������� ������������
double xmax;// ������� ������� ������������
int threads;// ���������� �������
double dmax;// ������������ ���������� ����������� ��������� �� ������ ����
double dm;// ������������ ��������� ����������� ��������� �� ������ ����
double v1;// �������� ������� ��������
double v2;// �������� ������� ��������
double x10;// ���������� ������ ������� ��������
double x20;// ���������� ������ ������� ��������
double xb; // ���������� ������-�������
double mu; // �������� ������-�������
double l; // �������� ������ ������-�������
double b;// �������� b (��������� �������)
double** f;// ������ �������� ����������
double* xmas;// ������ ����� �� ��� ������� � ����� � �������� �������� ������
double* ymas;// ������ ����� �� ��� ������� � ����� � �������� �������� ������
double xmin_kink;// ����������� �������� ��������� ������������
double xmax_kink;// ������������ �������� ��������� ������������
int size_kink;// ���������� ����� � ������� ������� �����
int divx;// �������� ���������� ����� �� ��� ��������� ��� ���������� � ����
int divt;// �������� ���������� ����� �� ��� ������� ��� ���������� � ����
string str;
real_1d_array xa;
real_1d_array ya;
barycentricinterpolant p;
spline1dinterpolant s;

// ��������������� ������� ��� ������� ����������� �����
void int_function_1_func(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(sqrt(1 + b*b*(1 - cos(x))) - 1);
}

// ��������, ������������ �������� ����� � ������� ����
double Fun(double alpha)
{
    autogkstate s;
    double v = 0.0;
    autogkreport rep;
    autogksmooth(pi, alpha, s);
    alglib::autogkintegrate(s, int_function_1_func);
    autogkresults(s, v, rep);
	return v;
}

// ��������� �������
// ��������� ����
double f0_single_kink(const double x, const double t, const double v1, const double x10, const double b)
{

	double f1 = 0.0;
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
	return f1;
}

// ������� ����
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

// ������� ������������ ����������
double delta_barrier(double x)
{
	return mu*l*exp(-sqr(l*(x - xb)))/sqrt(pi);
}

int main(int argc, char *argv[])
{
	//��������� ����� � �������
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
	//��������� ���������
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
	// ��������� ������-�������
	getline(in_file, str);
	getline(in_file, str);
	xb = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	mu = atof(str.data());
	in_file.close();
	l = 10;
	cout << "Begin computation." << endl;
	double begin = omp_get_wtime();
	// ���������� �������� ��������� �����
	// ������ ������������� ������� � ������� ������ �������� ����� �����
	vector<double> xmas_temp; 
	vector<double> ymas_temp;
	for(double alpha = 0.00001; alpha < 2*pi; alpha += 0.00005)
	{
		xmas_temp.push_back(Fun(alpha)/2);
		ymas_temp.push_back(alpha);
	}
	xmax_kink = max(xmas_temp[0], xmas_temp[xmas_temp.size()-1]);
	xmin_kink = min(xmas_temp[0], xmas_temp[xmas_temp.size()-1]);
	size_kink = xmas_temp.size();
	xmas = new double[size_kink];
	ymas = new double[size_kink];
	for(int i = 0; i < size_kink; i++)
	{
		xmas[i] = xmas_temp[i];
		ymas[i] = ymas_temp[i];
	}
	// ������ ������, �������������� �������� ��������� �����
	xa.setcontent(size_kink, xmas);
	ya.setcontent(size_kink, ymas);
	spline1dbuildlinear(xa, ya, s);

	//��������� ������� 
	hx = (xmax-xmin)/(nx-1);
	ht = hx/2;
	nt = (int)floor((tmax-tmin)/ht)+1; // ����� ����� ����� �� �������
	f = (double**)calloc(nx, sizeof(double));//����� ��������� �� ������� �������
	for(int i = 0; i < nx; i++)
	{
		f[i] = (double*)calloc(masnt, sizeof(double));//����� ��������� �� ������� �������
	};
	// ��������� ��������� �������
	// ���� ��������� ����
	if (v2 == 0.0)
		for (int x = 0; x < nx; x++)
		{
			f[x][0] = f0_single_kink(xmin+x*hx, tmin, v1, x10, b);
			f[x][1] = f0_single_kink(xmin+x*hx, tmin+ht, v1, x10, b);
		}
	else
		// ���� ���� �������
		for (int x = 0; x < nx; x++)
		{
			f[x][0] = f0_double_kink(xmin+x*hx, tmin, v1, v2, x10, x20, b);
			f[x][1] = f0_double_kink(xmin+x*hx, tmin+ht, v1, v2, x10, x20, b);
		};
	omp_set_num_threads(threads);
	// �������� �����
	#pragma omp parallel private (myid) shared (f, threads, masnt, nx, divt, divx, nt, mu, xb, l)
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
							- (1 + delta_barrier(xmin + x*hx))*b*b*ht*ht*sin(f[x][t-1])/sqrt(1 + b*b*(1 - cos(f[x][t-1])));
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
	cout << "Computation end in " << end - begin << " second(s)." << endl;
	cin.get();
	return 0;
}

