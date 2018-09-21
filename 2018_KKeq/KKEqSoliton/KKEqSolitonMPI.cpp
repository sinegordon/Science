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
int myid;// Номер процесса
int numprocs;// Число процессов
int now_nx;// Число точек данного процесса по пространственной оси
int nx;// Общее число точек по пространственной оси
int now_nt;// Число точек данного процесса по временной оси
double hx;// Шаг по пространственной оси
double ht;// Шаг по временной оси
double tmin;// Нижняя граница времени
double tmax;// Верхняя граница времени
double xmin;// Нижняя граница пространства
double xmax;// Верхняя граница пространства
double now_xmin;// Нижняя граница пространства для процесса
double eps;// Точность
double dmax;// Максимальная глобальная погрешность измерений на данном шаге
double dm;// Максимальная локальная погрешность измерений на данном шаге
double v1;// Скорость первого импульса
double v2;// Скорость второго импульса
double x10;// Координата центра первого начального импульса
double x20;// Координата центра второго начального импульса
double b;// Параметр b (отношение энергий)
double a; // Амплитуда ВЧ поля
double** f;// Массив значений потенциала
double* ff;// Массив значений потенциала для сохранения в файл (в процессе root)
double* fff;// Массив значений потенциала для передачи в процесс root
int* displs;// Сдвиг начала массивов-полос при пересылке и сборке
int* rcounts;// Количество пересылаемых элементов массивов-полос при пересылке и сборке
double* xmas;// Массив точек по оси абсцисс в файле с таблично заданным кинком
double* ymas;// Массив точек по оси ординат в файле с таблично заданным кинком
double* xmas1;// Массив точек по оси абсцисс в файле с таблично заданным током
double* ymas1;// Массив точек по оси ординат в файле с таблично заданным током
double xmin_kink;// Минимальное значение интервала интерполяции
double xmax_kink;// Максимальное значение интервала интерполяции
int size_kink;// Количество точек в таблице задания кинка (и тока)
int divx;// Делитель количества точек по оси координта при сохранении в файл
int divt;// Делитель количества точек по оси времени при сохранении в файл
double alpha_hf;// Вспомогательная переменная для потенциала при усреднении по периоду ВЧ поля
double F0; // F(0,a,b) в подкоренном выражении интеграла для усредненного кинка
string str;

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

//Переводим string в double
double StrToDouble(string str)
{
	return atof(str.data());
}

//Переводим string в int
int StrToInt(string str)
{
	return atoi(str.data());
}

// Подынтегральная функция при усредении по периоду ВЧ поля
void int_function_HF(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = sqrt(1 + b*b*(1 - cos(alpha_hf + a*sin(x))));
}

// Интеграл, определяющий значение функции в подкоренном выражении интеграла для кинка (F(alpha, a, b) в статье)
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

// Подынтегральная функция при неявном определении кинка усредненного по периоду ВЧ поля
void int_function_av_kink(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(F(x) - F0);
}

// Подынтегральная функция при неявном определении кинка
void int_function_1_func(double x, double xminusa, double bminusx, double &y, void *ptr) 
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
    alglib::autogkintegrate(s, int_function_1_func);
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

// Подынтегральная функция для тока усредненного по периоду ВЧ поля
void int_function_av_current(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
	y = b*b*sin(alpha_hf + a*sin(x))/sqrt(1 + b*b*(1 - cos(alpha_hf + a*sin(x))));
};

// Функция тока с усреднением по ВЧ полю
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

// Начальный профиль одиночного импульса
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
	// Инициализация параллельной библиотеки
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if (myid==0)
	{
		//Параметры сетки и решения
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
		// Амплитуда ВЧ поля
		getline(in_file, str);
		getline(in_file, str);
		a = atof(str.data());		
		in_file.close();
		cout << "Begin interpolation" << endl;
		// Построение таблично заданного кинка
		// Строим неравномерную таблицу с помощью неявно заданной формы кинка
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
		//Загрузка таблично заданного кинка
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
	
	// Рассылка начальных данных и параметров задачи по процессам
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
	// Строим сплайн, соответстующий таблично заданному кинку
	if(myid==0)
		cout << "Build splines" << endl;
	xa.setcontent(size_kink, xmas);
	ya.setcontent(size_kink, ymas);
	spline1dbuildlinear(xa, ya, s);
	xa1.setcontent(size_kink, xmas1);
	ya1.setcontent(size_kink, ymas1);
	spline1dbuildlinear(xa1, ya1, s1);

	//Начальные условия
	if(myid==0)
		cout << "Fill initial data" << endl;
	hx = (xmax-xmin)/(nx-1);
	ht = hx/2;
	if ((myid==0 || myid==numprocs-1) && numprocs > 1)
	{
		now_nx = (int)nx/numprocs+1; // Число узлов на полосу на крайних процессах
	};
	if (myid!=0&&myid!=(numprocs-1)&&numprocs>1)
	{
		now_nx = (int)nx/numprocs+2; // Число узлов на полосу на внутренних процессах
	};
	if (numprocs==1)
	{
		now_nx = nx; // Число узлов на полосу при одном процессе 
	};
	
	now_nt = (int)floor((tmax-tmin)/ht)+1; // Общее число узлов по времени
	f = (double**)calloc(now_nx, sizeof(double));//Число элементов по первому индексу
	for(int i = 0; i < now_nx; i++)
	{
		f[i] = (double*)calloc(now_nt, sizeof(double));//Число элементов по второму индексу
	};
	// Вычисляем минимальную границу по пространству для процесса
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
	// Заполняем начальные условия
	for (int x = 0; x < now_nx; x++)
	{
		f[x][0]=f0_single(now_xmin+x*hx, tmin, v1, x10);
		f[x][1]=f0_single(now_xmin+x*hx, tmin+ht, v1, x10);
	};

	// Включаем время
	// Первое приближение
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
						- ht*ht*spline1dcalc(s1, f[x][t-1]);//spline1dcalc - ток
		};
		f[0][t] = 2*f[1][t]-f[2][t];
		f[now_nx-1][t] = 2*f[now_nx-2][t]-f[now_nx-3][t];
	};
	//dmax=0;
	
	dmax=eps+1; // Максимальная невязка такая, чтобы начался цикл
	int tag = 1;
	double temp;
	MPI_Status status;
	if(numprocs > 1)
	{
		if(myid==0)
			cout << "Good approximation cycle is on" << endl;
		while(dmax > eps)
		{
			dm=0; // Текущая максимальная невязка
			// Посылаем следующему процессу предпоследнюю строку и принимаем от предыдущего
			// строку на первое место
			if(myid != numprocs - 1)
				MPI_Send(f[now_nx-2], now_nt, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD);		
			if(myid != 0)
				MPI_Recv(f[0], now_nt, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD, &status);
			// Посылаем предыдущему процессу вторую строку и принимаем от следующего
			// строку на первое место
			if(myid != 0)
				MPI_Send(f[1], now_nt, MPI_DOUBLE, myid-1, tag, MPI_COMM_WORLD);
			if(myid != numprocs-1)
				MPI_Recv(f[now_nx-1], now_nt, MPI_DOUBLE, myid+1, tag, MPI_COMM_WORLD, &status);
			// Следующая итерация
			for(int t=2;t<now_nt;t++)
			{
				for(int x=1;x<now_nx-1;x++)
				{
					temp = f[x][t];
					f[x][t] = (ht*ht)*(f[x-1][t-1]+f[x+1][t-1]-2*f[x][t-1])/(hx*hx)+2*f[x][t-1]-f[x][t-2]
							- ht*ht*spline1dcalc(s1, f[x][t-1]);//spline1dcalc - ток
					dm = fabs(temp-f[x][t]);
				};
			};		
			MPI_Reduce(&dm, &dmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			if (myid == 0)
				cout << "dmax = " << dmax << endl;
			MPI_Bcast(&dmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		};
	};
		
	// Переписываем вычисленный массив f[][] в одномерный fff[] и пересылаем его
	//в процесс root и там собираем в массив ff для сохранения в файл
	// Сборка
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
	// Пересылка
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
	// Вывод в файл потенциала
	if (myid==0)
	{
		cout << "End receving files" << endl;
		cout << "Saving files" << endl;
		out_file.open("out_f.sgo");
		/*out_file << 2 << endl;// Выводим сколько переменных
		out_file << "x" << endl;// Выводим название первой переменной
		out_file << "t" << endl;// Выводим название второй переменной
		out_file << nx << endl;// Выводим число точек по первой переменной (исключая нахлесты)
		out_file << now_nt << endl;// Выводим число точек по второй переменной
		out_file << hx << endl; // Выводим шаг по первой переменной
		out_file << ht << endl; // Выводим шаг по второй переменной
		out_file << xmin << endl; // Выводим нижнюю границу по первой переменной
		out_file << xmax << endl; // Выводим нижнюю границу по первой переменной
		out_file << tmin << endl; // Выводим нижнюю границу по второй переменной
		out_file << tmax << endl; // Выводим нижнюю границу по второй переменной
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
		// Вывод в файл напряженности
		/*
		out_file.open("out_e.sgo");
		out_file << 2 << endl;// Выводим сколько переменных
		out_file << "x" << endl;// Выводим название первой переменной
		out_file << "t" << endl;// Выводим название второй переменной
		out_file << nx-1 << endl;// Выводим число точек по первой переменной (исключая нахлесты)
		out_file << now_nt << endl;// Выводим число точек по второй переменной
		out_file << hx << endl; // Выводим шаг по первой переменной
		out_file << ht << endl; // Выводим шаг по второй переменной
		out_file << xmin << endl; // Выводим нижнюю границу по первой переменной
		out_file << xmax << endl; // Выводим нижнюю границу по первой переменной
		out_file << tmin << endl; // Выводим нижнюю границу по второй переменной
		out_file << tmax << endl; // Выводим нижнюю границу по второй переменной
		
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

