// Soliton.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include <stdlib.h>
#include <string>
#include <limits>
#include "alglib/src/interpolation.h"

#define pi 3.1415926535897932385

using namespace std;
using namespace alglib;

ifstream in_file;
ofstream out_file;
int myid;// ����� ��������
int numprocs;// ����� ���������
int now_nx;// ����� ����� ������� �������� �� ���������������� ���
int nx;// ����� ����� ����� �� ���������������� ���
int now_nt;// ����� ����� ������� �������� �� ��������� ���
double hx;// ��� �� ���������������� ���
double ht;// ��� �� ��������� ���
double tmin;// ������ ������� �������
double tmax;// ������� ������� �������
double xmin;// ������ ������� ������������
double xmax;// ������� ������� ������������
double now_xmin;// ������ ������� ������������ ��� ��������
double eps;// ��������
double dmax;// ������������ ���������� ����������� ��������� �� ������ ����
double dm;// ������������ ��������� ����������� ��������� �� ������ ����
double v1;// �������� ������� ��������
double v2;// �������� ������� ��������
double x10;// ���������� ������ ������� ���������� ��������
double x20;// ���������� ������ ������� ���������� ��������
double b;// �������� b (��������� �������)
double a; // ��������� �� ����
double** f;// ������ �������� ����������
double* ff;// ������ �������� ���������� ��� ���������� � ���� (� �������� root)
double* fff;// ������ �������� ���������� ��� �������� � ������� root
int* displs;// ����� ������ ��������-����� ��� ��������� � ������
int* rcounts;// ���������� ������������ ��������� ��������-����� ��� ��������� � ������
double* xmas;// ������ ����� �� ��� ������� � ����� � �������� �������� ������
double* ymas;// ������ ����� �� ��� ������� � ����� � �������� �������� ������
double* xmas1;// ������ ����� �� ��� ������� � ����� � �������� �������� �����
double* ymas1;// ������ ����� �� ��� ������� � ����� � �������� �������� �����
double xmin_kink;// ����������� �������� ��������� ������������
double xmax_kink;// ������������ �������� ��������� ������������
int size_kink;// ���������� ����� � ������� ������� ����� (� ����)
int divx;// �������� ���������� ����� �� ��� ��������� ��� ���������� � ����
int divt;// �������� ���������� ����� �� ��� ������� ��� ���������� � ����
double alpha_hf;// ��������������� ���������� ��� ���������� ��� ���������� �� ������� �� ����
double F0; // F(0,a,b) � ����������� ��������� ��������� ��� ������������ �����
string str;

// ��������� ��� ���������� ������� �����
real_1d_array xa;
real_1d_array ya;
barycentricinterpolant p;
spline1dinterpolant s;
// ��������� ��� ���������� ������� ����
real_1d_array xa1;
real_1d_array ya1;
barycentricinterpolant p1;
spline1dinterpolant s1;

//��������� string � double
double StrToDouble(string str)
{
	return atof(str.data());
}

//��������� string � int
int StrToInt(string str)
{
	return atoi(str.data());
}

// ��������������� ������� ��� ��������� �� ������� �� ����
void int_function_HF(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = sqrt(1 + b*b*(1 - cos(alpha_hf + a*sin(x))));
}

// ��������, ������������ �������� ������� � ����������� ��������� ��������� ��� ����� (F(alpha, a, b) � ������)
double F(double alpha)
{
    autogkstate s;
    double v = 0.0;
    autogkreport rep;
    autogksmooth(0, 2*pi, s);
	alpha_hf = alpha;
    alglib::autogkintegrate(s, int_function_HF);
    autogkresults(s, v, rep);
	return v/2.0/pi;
}

// ��������������� ������� ��� ������� ����������� ����� ������������ �� ������� �� ����
void int_function_av_kink(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(F(x) - F0);
}

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

// ��������, ������������ �������� ������������� ����� � ������� ����
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

// ��������������� ������� ��� ���� ������������ �� ������� �� ����
void int_function_av_current(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
	y = b*b*sin(alpha_hf + a*sin(x))/sqrt(1 + b*b*(1 - cos(alpha_hf + a*sin(x))));
};

// ������� ���� � ����������� �� �� ����
double currentHF(double f)
{
	autogkstate s;
	double v = 0.0;
	autogkreport rep;
	autogksmooth(0, 2*pi, s);
	alpha_hf = f;
	alglib::autogkintegrate(s, int_function_av_current);
	autogkresults(s, v, rep);
	return v/2.0/pi;
};

double current(double f)
{
	return b*b*sin(f)/sqrt(1 + b*b*(1 - cos(f)));
};

// ��������� ������� ���������� ��������
double f0_single(const double& x, const double& t, const double& v1, const double& x10)
{

	double f1 = 0.0;
	double arg1 = (x-x10-v1*t)/sqrt(1-v1*v1);

	if(arg1 < xmin_kink)
	{
		f1 = 0.0;
	}
	else
	{
		if (arg1 > xmax_kink)
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

int main(int argc, char *argv[])
{
	// ������������� ������������ ����������
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if (myid==0)
	{
		//��������� ����� � �������
		in_file.open("inMPI.txt");
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
		eps = atof(str.data());
		getline(in_file, str);
		getline(in_file, str);
		b = atof(str.data());
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
		// ��������� �� ����
		getline(in_file, str);
		getline(in_file, str);
		a = atof(str.data());		
		in_file.close();
		cout << "Begin interpolation" << endl;
		// ���������� �������� ��������� �����
		// ������ ������������� ������� � ������� ������ �������� ����� �����
		vector<double> xmas_temp; 
		vector<double> ymas_temp;
		F0 = F(0);
		int k = 0;
		for(double alpha = 0.000001; alpha < 2*pi; alpha += 0.0001)
		{
			if(k % 100 == 0)
				cout << "Compute alpha = " << alpha << endl;
			xmas_temp.push_back(FunHF(alpha)/2);
			ymas_temp.push_back(alpha);
			k += 1;
		}
		xmax_kink = max(xmas_temp[0], xmas_temp[xmas_temp.size()-1]);
		xmin_kink = min(xmas_temp[0], xmas_temp[xmas_temp.size()-1]);
		size_kink = xmas_temp.size();
		xmas = new double[size_kink];
		ymas = new double[size_kink];
		xmas1 = new double[size_kink];
		ymas1 = new double[size_kink];
		for(int i = 0; i < size_kink; i++)
		{
			xmas[i] = xmas_temp[i];
			ymas[i] = ymas_temp[i];
			xmas1[i] = ymas_temp[i];
			ymas1[i] = currentHF(xmas1[i]);
		}
		/*
		//�������� �������� ��������� �����
		in_file.open("kink.dat");
		getline(in_file, str);
		size_kink = atoi(str.data());
		xmin_kink = 0.0;
		xmax_kink = 0.0;
		xmas = new double[size_kink];
		ymas = new double[size_kink];
		for(int i = 0; i < size_kink; i++)
		{
			getline(in_file, str);
			xmas[i] = atof(str.data());
			if(xmas[i] < xmin_kink)
				xmin_kink = xmas[i];
			if(xmas[i] > xmax_kink)
				xmax_kink = xmas[i];
			getline(in_file, str);
			ymas[i] = atof(str.data());
		}
		in_file.close();
		*/
		cout << "Done interpolation" << endl;
		cout << "Begin solve wave equation" << endl;
	}
	
	// �������� ��������� ������ � ���������� ������ �� ���������
	if(myid==0)
		cout << "Broadcast data to processes" << endl;
 	MPI_Bcast(&nx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&xmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&v1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&v2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&x10, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&x20, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&xmin_kink, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&xmax_kink, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&size_kink, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myid > 0) 
	{
		xmas = new double[size_kink];
		ymas = new double[size_kink];
		xmas1 = new double[size_kink];
		ymas1 = new double[size_kink];
	}
	MPI_Bcast(xmas, size_kink, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ymas, size_kink, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(xmas1, size_kink, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ymas1, size_kink, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// ������ ������, �������������� �������� ��������� �����
	if(myid==0)
		cout << "Build splines" << endl;
	xa.setcontent(size_kink, xmas);
	ya.setcontent(size_kink, ymas);
	spline1dbuildlinear(xa, ya, s);
	xa1.setcontent(size_kink, xmas1);
	ya1.setcontent(size_kink, ymas1);
	spline1dbuildlinear(xa1, ya1, s1);

	//��������� �������
	if(myid==0)
		cout << "Fill initial data" << endl;
	hx = (xmax-xmin)/(nx-1);
	ht = hx/2;
	if ((myid==0 || myid==numprocs-1) && numprocs > 1)
	{
		now_nx = (int)nx/numprocs+1; // ����� ����� �� ������ �� ������� ���������
	};
	if (myid!=0&&myid!=(numprocs-1)&&numprocs>1)
	{
		now_nx = (int)nx/numprocs+2; // ����� ����� �� ������ �� ���������� ���������
	};
	if (numprocs==1)
	{
		now_nx = nx; // ����� ����� �� ������ ��� ����� �������� 
	};
	
	now_nt = (int)floor((tmax-tmin)/ht)+1; // ����� ����� ����� �� �������
	f = (double**)calloc(now_nx, sizeof(double));//����� ��������� �� ������� �������
	for(int i = 0; i < now_nx; i++)
	{
		f[i] = (double*)calloc(now_nt, sizeof(double));//����� ��������� �� ������� �������
	};
	// ��������� ����������� ������� �� ������������ ��� ��������
	if (myid==0)
	{
		now_xmin = xmin;
	};
	if (myid==numprocs-1)
	{
		now_xmin = xmax-(now_nx-1)*hx;
	};
	if ((myid!=0)&&(myid!=numprocs-1))
	{
		now_xmin = xmin + nx/numprocs*myid*hx;
	};
	// ��������� ��������� �������
	for (int x = 0; x < now_nx; x++)
	{
		f[x][0]=f0_single(now_xmin+x*hx, tmin, v1, x10);
		f[x][1]=f0_single(now_xmin+x*hx, tmin+ht, v1, x10);
	};

	// �������� �����
	// ������ �����������
	if(myid==0)
		cout << "Time is on" << endl;
	if(myid==0)
		cout << "First approximation" << endl;
	for(int t = 2; t < now_nt; t++)
	{
		if(myid==0 && t%100==0)
			cout << "Time = " << t << " from " << now_nt << endl;
		for(int x = 1; x < now_nx-1; x++)
		{
			//f[x][t]=f0(now_xmin+x*hx, tmin+t*ht, v1, v2, x10, x20, b);
			f[x][t] = (ht*ht)*(f[x-1][t-1]+f[x+1][t-1]-2*f[x][t-1])/(hx*hx)+2*f[x][t-1]-f[x][t-2]
						- ht*ht*spline1dcalc(s1, f[x][t-1]);//spline1dcalc - ���
		};
		f[0][t] = 2*f[1][t]-f[2][t];
		f[now_nx-1][t] = 2*f[now_nx-2][t]-f[now_nx-3][t];
	};
	//dmax=0;
	
	dmax=eps+1; // ������������ ������� �����, ����� ������� ����
	int tag = 1;
	double temp;
	MPI_Status status;
	if(numprocs > 1)
	{
		if(myid==0)
			cout << "Good approximation cycle is on" << endl;
		while(dmax > eps)
		{
			dm=0; // ������� ������������ �������
			// �������� ���������� �������� ������������� ������ � ��������� �� �����������
			// ������ �� ������ �����
			if(myid != numprocs - 1)
				MPI_Send(f[now_nx-2], now_nt, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);		
			if(myid != 0)
				MPI_Recv(f[0], now_nt, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
			// �������� ����������� �������� ������ ������ � ��������� �� ����������
			// ������ �� ������ �����
			if(myid != 0)
				MPI_Send(f[1], now_nt, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
			if(myid != numprocs-1)
				MPI_Recv(f[now_nx-1], now_nt, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
			// ��������� ��������
			for(int t=2;t<now_nt;t++)
			{
				for(int x=1;x<now_nx-1;x++)
				{
					temp = f[x][t];
					f[x][t] = (ht*ht)*(f[x-1][t-1]+f[x+1][t-1]-2*f[x][t-1])/(hx*hx)+2*f[x][t-1]-f[x][t-2]
							- ht*ht*spline1dcalc(s1, f[x][t-1]);//spline1dcalc - ���
					dm = fabs(temp-f[x][t]);
				};
			};		
			MPI_Reduce(&dm, &dmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			if (myid == 0)
				cout << "dmax = " << dmax << endl;
			MPI_Bcast(&dmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		};
	};
		
	// ������������ ����������� ������ f[][] � ���������� fff[] � ���������� ���
	//� ������� root � ��� �������� � ������ ff ��� ���������� � ����
	// ������
	if (myid==0)
	{
		cout << "Done solve wave equation" << endl;
		cout << "Receving files" << endl;
		ff = (double*)calloc(nx*now_nt, sizeof(double));
		fff = (double*)calloc(now_nx*now_nt, sizeof(double));
		for(int x = 0; x < now_nx; x++)
		{
			for(int t = 0; t < now_nt; t++)
			{
				fff[x*now_nt+t] = f[x][t];
			};
		};
		displs = (int*)calloc(numprocs, sizeof(int));
		rcounts = (int*)calloc(numprocs, sizeof(int));
		for (int i = 0; i < numprocs; i++) 
		{
			if (i==0)
			{
				displs[i] = 0;
		        rcounts[i] = now_nx*now_nt;
			}
			if (i==numprocs-1)
			{
				displs[i] = displs[i-1]+rcounts[i-1];
				rcounts[i] = (now_nx-2)*now_nt;
				//ff = &ff[2];
			};
			if ((i!=numprocs-1)&&(i!=0))
			{
				displs[i] = displs[i-1]+rcounts[i-1];
				rcounts[i] = (now_nx-1)*now_nt;
				//ff = &ff[1];
			};
		};
	}
	else
	{
		fff = (double*)calloc((now_nx-2)*now_nt, sizeof(double));
		for(int x = 0; x < now_nx-2; x++)
		{
			for(int t = 0; t < now_nt; t++)
			{
				fff[x*now_nt+t] = f[x+2][t];
			};
		};
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// ���������
	if (numprocs==1)
	{
		MPI_Gather(fff, now_nx*now_nt, MPI_DOUBLE, ff, now_nx*now_nt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else
	{
		if(myid==0)
		{
			MPI_Gatherv(fff, now_nx*now_nt, MPI_DOUBLE, ff, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Gatherv(fff, (now_nx-2)*now_nt, MPI_DOUBLE, ff, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		};
	};
	MPI_Finalize();
	// ����� � ���� ����������
	if (myid==0)
	{
		cout << "End receving files" << endl;
		cout << "Saving files" << endl;
		out_file.open("out_f.sgo");
		/*out_file << 2 << endl;// ������� ������� ����������
		out_file << "x" << endl;// ������� �������� ������ ����������
		out_file << "t" << endl;// ������� �������� ������ ����������
		out_file << nx << endl;// ������� ����� ����� �� ������ ���������� (�������� ��������)
		out_file << now_nt << endl;// ������� ����� ����� �� ������ ����������
		out_file << hx << endl; // ������� ��� �� ������ ����������
		out_file << ht << endl; // ������� ��� �� ������ ����������
		out_file << xmin << endl; // ������� ������ ������� �� ������ ����������
		out_file << xmax << endl; // ������� ������ ������� �� ������ ����������
		out_file << tmin << endl; // ������� ������ ������� �� ������ ����������
		out_file << tmax << endl; // ������� ������ ������� �� ������ ����������
		*/
		for(int t = 0; t < now_nt; t+=divt)
		{
			for(int x = 0; x < nx - 1; x+=divx)
			{
				out_file << ff[x*now_nt+t] << " ";//ff[x*now_nt+t] << endl;
			};
			out_file << ff[(nx - 1)*now_nt+t] << endl;
		};
		out_file.close();
		// ����� � ���� �������������
		/*
		out_file.open("out_e.sgo");
		out_file << 2 << endl;// ������� ������� ����������
		out_file << "x" << endl;// ������� �������� ������ ����������
		out_file << "t" << endl;// ������� �������� ������ ����������
		out_file << nx-1 << endl;// ������� ����� ����� �� ������ ���������� (�������� ��������)
		out_file << now_nt << endl;// ������� ����� ����� �� ������ ����������
		out_file << hx << endl; // ������� ��� �� ������ ����������
		out_file << ht << endl; // ������� ��� �� ������ ����������
		out_file << xmin << endl; // ������� ������ ������� �� ������ ����������
		out_file << xmax << endl; // ������� ������ ������� �� ������ ����������
		out_file << tmin << endl; // ������� ������ ������� �� ������ ����������
		out_file << tmax << endl; // ������� ������ ������� �� ������ ����������
		
		for(int x=0;x<nx-1;x++)
		{
			for(int t=0;t<now_nt;t++)
			{
				out_file << (ff[(x+1)*now_nt+t]-ff[x*now_nt+t])/hx << endl;
			};
		};
		out_file.close();
		for(int t = 0; t < now_nt; t++)
		{
			for(int x = 0; x < nx - 1; x++)
			{
				out_file << (ff[(x+1)*now_nt+t]-ff[x*now_nt+t])/hx << " ";//ff[x*now_nt+t] << endl;
			};
			out_file << (ff[(nx - 2)*now_nt+t]-ff[(nx - 1)*now_nt+t])/hx << endl;
		};
		out_file.close();
		*/
		cout << "End saving files" << endl;
	};
	return 0;
}

