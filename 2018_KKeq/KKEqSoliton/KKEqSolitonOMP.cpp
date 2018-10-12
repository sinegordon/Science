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
<<<<<<< HEAD
int myid;// –ù–æ–º–µ—Ä –ø–æ—Ç–æ–∫–∞
int numthreads;// –ß–∏—Å–ª–æ –ø–æ—Ç–æ–∫–æ–≤
int nx;// –û–±—â–µ–µ —á–∏—Å–ª–æ —Ç–æ—á–µ–∫ –ø–æ –ø—Ä–æ—Å—Ç—Ä–∞–Ω—Å—Ç–≤–µ–Ω–Ω–æ–π –æ—Å–∏
int nt;// –û–±—â–µ–µ —á–∏—Å–ª–æ —Ç–æ—á–µ–∫ –ø–æ –≤—Ä–µ–º–µ–Ω–Ω–æ–π –æ—Å–∏
int masnt;// –ß–∏—Å–ª–æ —Ç–æ—á–µ–∫ –≤ –º–∞—Å—Å–∏–≤–µ –ø–æ –≤—Ä–µ–º–µ–Ω–Ω–æ–π –æ—Å–∏
double hx;// –®–∞–≥ –ø–æ –ø—Ä–æ—Å—Ç—Ä–∞–Ω—Å—Ç–≤–µ–Ω–Ω–æ–π –æ—Å–∏
double ht;// –®–∞–≥ –ø–æ –≤—Ä–µ–º–µ–Ω–Ω–æ–π –æ—Å–∏
double tmin;// –ù–∏–∂–Ω—è—è –≥—Ä–∞–Ω–∏—Ü–∞ –≤—Ä–µ–º–µ–Ω–∏
double tmax;// –í–µ—Ä—Ö–Ω—è—è –≥—Ä–∞–Ω–∏—Ü–∞ –≤—Ä–µ–º–µ–Ω–∏
double xmin;// –ù–∏–∂–Ω—è—è –≥—Ä–∞–Ω–∏—Ü–∞ –ø—Ä–æ—Å—Ç—Ä–∞–Ω—Å—Ç–≤–∞
double xmax;// –í–µ—Ä—Ö–Ω—è—è –≥—Ä–∞–Ω–∏—Ü–∞ –ø—Ä–æ—Å—Ç—Ä–∞–Ω—Å—Ç–≤–∞
int threads;// –ö–æ–ª–∏—á–µ—Å—Ç–≤–æ –ø–æ—Ç–æ–∫–æ–≤
double dmax;// –ú–∞–∫—Å–∏–º–∞–ª—å–Ω–∞—è –≥–ª–æ–±–∞–ª—å–Ω–∞—è –ø–æ–≥—Ä–µ—à–Ω–æ—Å—Ç—å –∏–∑–º–µ—Ä–µ–Ω–∏–π –Ω–∞ –¥–∞–Ω–Ω–æ–º —à–∞–≥–µ
double dm;// –ú–∞–∫—Å–∏–º–∞–ª—å–Ω–∞—è –ª–æ–∫–∞–ª—å–Ω–∞—è –ø–æ–≥—Ä–µ—à–Ω–æ—Å—Ç—å –∏–∑–º–µ—Ä–µ–Ω–∏–π –Ω–∞ –¥–∞–Ω–Ω–æ–º —à–∞–≥–µ
double v1;// –°–∫–æ—Ä–æ—Å—Ç—å –ø–µ—Ä–≤–æ–≥–æ –∏–º–ø—É–ª—å—Å–∞
double v2;// –°–∫–æ—Ä–æ—Å—Ç—å –≤—Ç–æ—Ä–æ–≥–æ –∏–º–ø—É–ª—å—Å–∞
double x10;// –ö–æ–æ—Ä–¥–∏–Ω–∞—Ç–∞ —Ü–µ–Ω—Ç—Ä–∞ –ø–µ—Ä–≤–æ–≥–æ –∏–º–ø—É–ª—å—Å–∞
double x20;// –ö–æ–æ—Ä–¥–∏–Ω–∞—Ç–∞ —Ü–µ–Ω—Ç—Ä–∞ –≤—Ç–æ—Ä–æ–≥–æ –∏–º–ø—É–ª—å—Å–∞
double xb; // –ö–æ–æ—Ä–¥–∏–Ω–∞—Ç–∞ –¥–µ–ª—å—Ç–∞-–±–∞—Ä—å–µ—Ä–∞
double mu; // –ú–æ—â–Ω–æ—Å—Ç—å –¥–µ–ª—å—Ç–∞-–±–∞—Ä—å–µ—Ä–∞
double l; // –û–±—Ä–∞—Ç–Ω–∞—è —à–∏—Ä–∏–Ω–∞ –¥–µ–ª—å—Ç–∞-–±–∞—Ä—å–µ—Ä–∞
double a;// –ê–º–ø–ª–∏—Ç—É–¥–∞ –í–ß-–ø–æ–ª—è
double b;// –ü–∞—Ä–∞–º–µ—Ç—Ä b (–æ—Ç–Ω–æ—à–µ–Ω–∏–µ —ç–Ω–µ—Ä–≥–∏–π)
double** f;// –ú–∞—Å—Å–∏–≤ –∑–Ω–∞—á–µ–Ω–∏–π –ø–æ—Ç–µ–Ω—Ü–∏–∞–ª–∞
double* xmas;// –ú–∞—Å—Å–∏–≤ —Ç–æ—á–µ–∫ –ø–æ –æ—Å–∏ –∞–±—Å—Ü–∏—Å—Å –≤ —Ñ–∞–π–ª–µ —Å —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω—ã–º –∫–∏–Ω–∫–æ–º
double* ymas;// –ú–∞—Å—Å–∏–≤ —Ç–æ—á–µ–∫ –ø–æ –æ—Å–∏ –æ—Ä–¥–∏–Ω–∞—Ç –≤ —Ñ–∞–π–ª–µ —Å —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω—ã–º –∫–∏–Ω–∫–æ–º
double xmin_kink;// –ú–∏–Ω–∏–º–∞–ª—å–Ω–æ–µ –∑–Ω–∞—á–µ–Ω–∏–µ –∏–Ω—Ç–µ—Ä–≤–∞–ª–∞ –∏–Ω—Ç–µ—Ä–ø–æ–ª—è—Ü–∏–∏
double xmax_kink;// –ú–∞–∫—Å–∏–º–∞–ª—å–Ω–æ–µ –∑–Ω–∞—á–µ–Ω–∏–µ –∏–Ω—Ç–µ—Ä–≤–∞–ª–∞ –∏–Ω—Ç–µ—Ä–ø–æ–ª—è—Ü–∏–∏
double* xmas1;// –ú–∞—Å—Å–∏–≤ —Ç–æ—á–µ–∫ –ø–æ –æ—Å–∏ –∞–±—Å—Ü–∏—Å—Å –≤ —Ñ–∞–π–ª–µ —Å —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω—ã–º —Ç–æ–∫–æ–º
double* ymas1;// –ú–∞—Å—Å–∏–≤ —Ç–æ—á–µ–∫ –ø–æ –æ—Å–∏ –æ—Ä–¥–∏–Ω–∞—Ç –≤ —Ñ–∞–π–ª–µ —Å —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω—ã–º —Ç–æ–∫–æ–º
double fmin_kink;// –ú–∏–Ω–∏–º–∞–ª—å–Ω–æ–µ –∑–Ω–∞—á–µ–Ω–∏–µ –ø–æ—Ç–µ–Ω—Ü–∏–∞–ª–∞ –ø—Ä–∏ –∏–Ω—Ç–µ—Ä–ø–æ–ª—è—Ü–∏–∏
double fmax_kink;// –ú–∞–∫—Å–∏–º–∞–ª—å–Ω–æ–µ –∑–Ω–∞—á–µ–Ω–∏–µ –ø–æ—Ç–µ–Ω—Ü–∏–∞–ª–∞ –ø—Ä–∏ –∏–Ω—Ç–µ—Ä–ø–æ–ª—è—Ü–∏–∏
int size_kink;// –ö–æ–ª–∏—á–µ—Å—Ç–≤–æ —Ç–æ—á–µ–∫ –≤ —Ç–∞–±–ª–∏—Ü–µ –∑–∞–¥–∞–Ω–∏—è –∫–∏–Ω–∫–∞
int divx;// –î–µ–ª–∏—Ç–µ–ª—å –∫–æ–ª–∏—á–µ—Å—Ç–≤–∞ —Ç–æ—á–µ–∫ –ø–æ –æ—Å–∏ –∫–æ–æ—Ä–¥–∏–Ω–∞—Ç–∞ –ø—Ä–∏ —Å–æ—Ö—Ä–∞–Ω–µ–Ω–∏–∏ –≤ —Ñ–∞–π–ª
int divt;// –î–µ–ª–∏—Ç–µ–ª—å –∫–æ–ª–∏—á–µ—Å—Ç–≤–∞ —Ç–æ—á–µ–∫ –ø–æ –æ—Å–∏ –≤—Ä–µ–º–µ–Ω–∏ –ø—Ä–∏ —Å–æ—Ö—Ä–∞–Ω–µ–Ω–∏–∏ –≤ —Ñ–∞–π–ª
double betta;// –ö–æ—ç—Ñ—Ñ–∏—Ü–∏–µ–Ω—Ç —Ç—Ä–µ–Ω–∏—è
=======
int myid;// ÕÓÏÂ ÔÓÚÓÍ‡
int numthreads;// ◊ËÒÎÓ ÔÓÚÓÍÓ‚
int nx;// Œ·˘ÂÂ ˜ËÒÎÓ ÚÓ˜ÂÍ ÔÓ ÔÓÒÚ‡ÌÒÚ‚ÂÌÌÓÈ ÓÒË
int nt;// Œ·˘ÂÂ ˜ËÒÎÓ ÚÓ˜ÂÍ ÔÓ ‚ÂÏÂÌÌÓÈ ÓÒË
int masnt;// ◊ËÒÎÓ ÚÓ˜ÂÍ ‚ Ï‡ÒÒË‚Â ÔÓ ‚ÂÏÂÌÌÓÈ ÓÒË
double hx;// ÿ‡„ ÔÓ ÔÓÒÚ‡ÌÒÚ‚ÂÌÌÓÈ ÓÒË
double ht;// ÿ‡„ ÔÓ ‚ÂÏÂÌÌÓÈ ÓÒË
double tmin;// ÕËÊÌˇˇ „‡ÌËˆ‡ ‚ÂÏÂÌË
double tmax;// ¬ÂıÌˇˇ „‡ÌËˆ‡ ‚ÂÏÂÌË
double xmin;// ÕËÊÌˇˇ „‡ÌËˆ‡ ÔÓÒÚ‡ÌÒÚ‚‡
double xmax;// ¬ÂıÌˇˇ „‡ÌËˆ‡ ÔÓÒÚ‡ÌÒÚ‚‡
int threads;//  ÓÎË˜ÂÒÚ‚Ó ÔÓÚÓÍÓ‚
double dmax;// Ã‡ÍÒËÏ‡Î¸Ì‡ˇ „ÎÓ·‡Î¸Ì‡ˇ ÔÓ„Â¯ÌÓÒÚ¸ ËÁÏÂÂÌËÈ Ì‡ ‰‡ÌÌÓÏ ¯‡„Â
double dm;// Ã‡ÍÒËÏ‡Î¸Ì‡ˇ ÎÓÍ‡Î¸Ì‡ˇ ÔÓ„Â¯ÌÓÒÚ¸ ËÁÏÂÂÌËÈ Ì‡ ‰‡ÌÌÓÏ ¯‡„Â
double v1;// —ÍÓÓÒÚ¸ ÔÂ‚Ó„Ó ËÏÔÛÎ¸Ò‡
double v2;// —ÍÓÓÒÚ¸ ‚ÚÓÓ„Ó ËÏÔÛÎ¸Ò‡
double x10;//  ÓÓ‰ËÌ‡Ú‡ ˆÂÌÚ‡ ÔÂ‚Ó„Ó ËÏÔÛÎ¸Ò‡
double x20;//  ÓÓ‰ËÌ‡Ú‡ ˆÂÌÚ‡ ‚ÚÓÓ„Ó ËÏÔÛÎ¸Ò‡
double xb; //  ÓÓ‰ËÌ‡Ú‡ ‰ÂÎ¸Ú‡-·‡¸Â‡
double mu; // ÃÓ˘ÌÓÒÚ¸ ‰ÂÎ¸Ú‡-·‡¸Â‡
double l; // Œ·‡ÚÌ‡ˇ ¯ËËÌ‡ ‰ÂÎ¸Ú‡-·‡¸Â‡
double a;// ¿ÏÔÎËÚÛ‰‡ ¬◊-ÔÓÎˇ
double b;// œ‡‡ÏÂÚ b (ÓÚÌÓ¯ÂÌËÂ ˝ÌÂ„ËÈ)
double** f;// Ã‡ÒÒË‚ ÁÌ‡˜ÂÌËÈ ÔÓÚÂÌˆË‡Î‡
double* xmas;// Ã‡ÒÒË‚ ÚÓ˜ÂÍ ÔÓ ÓÒË ‡·ÒˆËÒÒ ‚ Ù‡ÈÎÂ Ò Ú‡·ÎË˜ÌÓ Á‡‰‡ÌÌ˚Ï ÍËÌÍÓÏ
double* ymas;// Ã‡ÒÒË‚ ÚÓ˜ÂÍ ÔÓ ÓÒË Ó‰ËÌ‡Ú ‚ Ù‡ÈÎÂ Ò Ú‡·ÎË˜ÌÓ Á‡‰‡ÌÌ˚Ï ÍËÌÍÓÏ
double xmin_kink;// ÃËÌËÏ‡Î¸ÌÓÂ ÁÌ‡˜ÂÌËÂ ËÌÚÂ‚‡Î‡ ËÌÚÂÔÓÎˇˆËË
double xmax_kink;// Ã‡ÍÒËÏ‡Î¸ÌÓÂ ÁÌ‡˜ÂÌËÂ ËÌÚÂ‚‡Î‡ ËÌÚÂÔÓÎˇˆËË
double* xmas1;// Ã‡ÒÒË‚ ÚÓ˜ÂÍ ÔÓ ÓÒË ‡·ÒˆËÒÒ ‚ Ù‡ÈÎÂ Ò Ú‡·ÎË˜ÌÓ Á‡‰‡ÌÌ˚Ï ÚÓÍÓÏ
double* ymas1;// Ã‡ÒÒË‚ ÚÓ˜ÂÍ ÔÓ ÓÒË Ó‰ËÌ‡Ú ‚ Ù‡ÈÎÂ Ò Ú‡·ÎË˜ÌÓ Á‡‰‡ÌÌ˚Ï ÚÓÍÓÏ
double fmin_kink;// ÃËÌËÏ‡Î¸ÌÓÂ ÁÌ‡˜ÂÌËÂ ÔÓÚÂÌˆË‡Î‡ ÔË ËÌÚÂÔÓÎˇˆËË
double fmax_kink;// Ã‡ÍÒËÏ‡Î¸ÌÓÂ ÁÌ‡˜ÂÌËÂ ÔÓÚÂÌˆË‡Î‡ ÔË ËÌÚÂÔÓÎˇˆËË
int size_kink;//  ÓÎË˜ÂÒÚ‚Ó ÚÓ˜ÂÍ ‚ Ú‡·ÎËˆÂ Á‡‰‡ÌËˇ ÍËÌÍ‡
int divx;// ƒÂÎËÚÂÎ¸ ÍÓÎË˜ÂÒÚ‚‡ ÚÓ˜ÂÍ ÔÓ ÓÒË ÍÓÓ‰ËÌ‡Ú‡ ÔË ÒÓı‡ÌÂÌËË ‚ Ù‡ÈÎ
int divt;// ƒÂÎËÚÂÎ¸ ÍÓÎË˜ÂÒÚ‚‡ ÚÓ˜ÂÍ ÔÓ ÓÒË ‚ÂÏÂÌË ÔË ÒÓı‡ÌÂÌËË ‚ Ù‡ÈÎ
double betta;//  Ó˝ÙÙËˆËÂÌÚ ÚÂÌËˇ
double w;// ◊‡ÒÚÓÚ‡ ÔÓÎˇ
double tau;// ’‡‡ÍÚÂÌÓÂ ‚ÂÏˇ ‚ÍÎ˛˜ÂÌËˇ ÔÓÎˇ
>>>>>>> 3ab42d3d87c37acc9d9821708d293f81c0f301c4
string str;
double * alpha_hf;// –í—Å–ø–æ–º–æ–≥–∞—Ç–µ–ª—å–Ω–∞—è –ø–µ—Ä–µ–º–µ–Ω–Ω–∞—è –¥–ª—è –ø–æ—Ç–µ–Ω—Ü–∏–∞–ª–∞ –ø—Ä–∏ —É—Å—Ä–µ–¥–Ω–µ–Ω–∏–∏ –ø–æ –ø–µ—Ä–∏–æ–¥—É –í–ß –ø–æ–ª—è
double F0; //  F(0,a,b) –≤ –ø–æ–¥–∫–æ—Ä–µ–Ω–Ω–æ–º –≤—ã—Ä–∞–∂–µ–Ω–∏–∏ –∏–Ω—Ç–µ–≥—Ä–∞–ª–∞ –¥–ª—è —É—Å—Ä–µ–¥–Ω–µ–Ω–Ω–æ–≥–æ –∫–∏–Ω–∫–∞
// –ü—Ä–µ–º–µ–Ω–Ω—ã–µ –¥–ª—è –ø–æ—Å—Ç—Ä–æ–µ–Ω–∏—è —Å–ø–ª–∞–π–Ω–∞ –∫–∏–Ω–∫–∞
real_1d_array xa;
real_1d_array ya;
barycentricinterpolant p;
spline1dinterpolant s;
// –ü—Ä–µ–º–µ–Ω–Ω—ã–µ –¥–ª—è –ø–æ—Å—Ç—Ä–æ–µ–Ω–∏—è —Å–ø–ª–∞–π–Ω–∞ —Ç–æ–∫–∞
real_1d_array xa1;
real_1d_array ya1;
barycentricinterpolant p1;
spline1dinterpolant s1;

// –ü–æ–¥—ã–Ω—Ç–µ–≥—Ä–∞–ª—å–Ω–∞—è —Ñ—É–Ω–∫—Ü–∏—è –ø—Ä–∏ —É—Å—Ä–µ–¥–µ–Ω–∏–∏ –ø–æ –ø–µ—Ä–∏–æ–¥—É –í–ß –ø–æ–ª—è
void int_function_HF(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = sqrt(1 + b*b*(1 - cos(alpha_hf[myid] + a*sin(x))));
}

// –ò–Ω—Ç–µ–≥—Ä–∞–ª, –æ–ø—Ä–µ–¥–µ–ª—è—é—â–∏–π –∑–Ω–∞—á–µ–Ω–∏–µ —Ñ—É–Ω–∫—Ü–∏–∏ –≤ –ø–æ–¥–∫–æ—Ä–µ–Ω–Ω–æ–º –≤—ã—Ä–∞–∂–µ–Ω–∏–∏ –∏–Ω—Ç–µ–≥—Ä–∞–ª–∞ –¥–ª—è –∫–∏–Ω–∫–∞ (F(alpha, a, b) –≤ —Å—Ç–∞—Ç—å–µ)
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

// –ü–æ–¥—ã–Ω—Ç–µ–≥—Ä–∞–ª—å–Ω–∞—è —Ñ—É–Ω–∫—Ü–∏—è –ø—Ä–∏ –Ω–µ—è–≤–Ω–æ–º –æ–ø—Ä–µ–¥–µ–ª–µ–Ω–∏–∏ –∫–∏–Ω–∫–∞ —É—Å—Ä–µ–¥–Ω–µ–Ω–Ω–æ–≥–æ –ø–æ –ø–µ—Ä–∏–æ–¥—É –í–ß –ø–æ–ª—è
void int_function_av_kink(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(F(x) - F0);
}

// –ü–æ–¥—ã–Ω—Ç–µ–≥—Ä–∞–ª—å–Ω–∞—è —Ñ—É–Ω–∫—Ü–∏—è –ø—Ä–∏ –Ω–µ—è–≤–Ω–æ–º –æ–ø—Ä–µ–¥–µ–ª–µ–Ω–∏–∏ –∫–∏–Ω–∫–∞
void int_function_kink(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
    y = 1/sqrt(sqrt(1 + b*b*(1 - cos(x))) - 1);
}

// –ò–Ω—Ç–µ–≥—Ä–∞–ª, –æ–ø—Ä–µ–¥–µ–ª—è—é—â–∏–π –∑–Ω–∞—á–µ–Ω–∏–µ –∫–∏–Ω–∫–∞ –≤ –Ω–µ—è–≤–Ω–æ–º –≤–∏–¥–µ
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

// –ò–Ω—Ç–µ–≥—Ä–∞–ª, –æ–ø—Ä–µ–¥–µ–ª—è—é—â–∏–π –∑–Ω–∞—á–µ–Ω–∏–µ —É—Å—Ä–µ–¥–µ–Ω–µ–Ω–Ω–æ–≥–æ –∫–∏–Ω–∫–∞ –≤ –Ω–µ—è–≤–Ω–æ–º –≤–∏–¥–µ
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

// –ù–∞—á–∞–ª—å–Ω—ã–π –ø—Ä–æ—Ñ–∏–ª—å
// –û–¥–∏–Ω–æ—á–Ω—ã–π –∫–∏–Ω–∫
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

// –î–≤–æ–π–Ω–æ–π –∫–∏–Ω–∫
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

// –§—É–Ω–∫—Ü–∏—è —Å—Ç—Ä—É–∫—Ç—É—Ä–Ω–æ–≥–æ –≤–æ–∑–º—É—â–µ–Ω–∏—è
inline double delta_barrier(double x)
{
	return mu*l*exp(-sqr(l*(x - xb)))/sqrt(pi);
}


// –ü–æ–¥—ã–Ω—Ç–µ–≥—Ä–∞–ª—å–Ω–∞—è —Ñ—É–Ω–∫—Ü–∏—è –¥–ª—è —Ç–æ–∫–∞ —É—Å—Ä–µ–¥–Ω–µ–Ω–Ω–æ–≥–æ –ø–æ –ø–µ—Ä–∏–æ–¥—É –í–ß –ø–æ–ª—è
void int_function_av_current(double x, double xminusa, double bminusx, double &y, void *ptr) 
{
	y = b*b*sin(alpha_hf[myid] + a*sin(x))/sqrt(1 + b*b*(1 - cos(alpha_hf[myid] + a*sin(x))));
};

// –§—É–Ω–∫—Ü–∏—è —Ç–æ–∫–∞ —Å —É—Å—Ä–µ–¥–Ω–µ–Ω–∏–µ–º –ø–æ –í–ß –ø–æ–ª—é
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
	//–ü–∞—Ä–∞–º–µ—Ç—Ä—ã —Å–µ—Ç–∫–∏ –∏ —Ä–µ—à–µ–Ω–∏—è
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
	//–ü–∞—Ä–∞–º–µ—Ç—Ä—ã –∏–º–ø—É–ª—å—Å–æ–≤
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
	// –ü–∞—Ä–∞–º–µ—Ç—Ä—ã –¥–µ–ª—å—Ç–∞-–±–∞—Ä—å–µ—Ä–∞
	getline(in_file, str);
	getline(in_file, str);
	xb = atof(str.data());	
	getline(in_file, str);
	getline(in_file, str);
	mu = atof(str.data());
	// –ê–º–ø–ª–∏—Ç—É–¥–∞ –í–ß –ø–æ–ª—è
	getline(in_file, str);
	getline(in_file, str);
	a = atof(str.data());
	// –ö–æ—ç—Ñ—Ñ–∏—Ü–∏–µ–Ω—Ç —Ç—Ä–µ–Ω–∏—è
	getline(in_file, str);
	getline(in_file, str);
	betta = atof(str.data());
	// ◊‡ÒÚÓÚ‡ ÔÓÎˇ
	getline(in_file, str);
	getline(in_file, str);
	w = atof(str.data());
	// ◊‡ÒÚÓÚ‡ ÔÓÎˇ
	getline(in_file, str);
	getline(in_file, str);
	tau = atof(str.data());

	in_file.close();
	l = 10;
	cout << "Begin interpolation kink and current" << endl;
	alpha_hf = new double[threads];
	myid = 0;
	F0 = F(0);
	double begin = omp_get_wtime();
	// –ü–æ—Å—Ç—Ä–æ–µ–Ω–∏–µ —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω–æ–≥–æ –∫–∏–Ω–∫–∞
	// –°—Ç—Ä–æ–∏–º –Ω–µ—Ä–∞–≤–Ω–æ–º–µ—Ä–Ω—É—é —Ç–∞–±–ª–∏—Ü—É —Å –ø–æ–º–æ—â—å—é –Ω–µ—è–≤–Ω–æ –∑–∞–¥–∞–Ω–Ω–æ–π —Ñ–æ—Ä–º—ã –∫–∏–Ω–∫–∞
	vector<double> xmas_temp; 
	vector<double> ymas_temp;
	int k = 0;
	for(double alpha = 0.000001; alpha < 2*pi; alpha += 0.0005)
	{
		xmas_temp.push_back(Fun(alpha)/2);
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
	// –°—Ç—Ä–æ–∏–º —Å–ø–ª–∞–π–Ω, —Å–æ–æ—Ç–≤–µ—Ç—Å—Ç–≤—É—é—â–∏–π —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω–æ–º—É –∫–∏–Ω–∫—É
	xa.setcontent(size_kink, xmas);
	ya.setcontent(size_kink, ymas);
	spline1dbuildcubic(xa, ya, s);
	// –°—Ç—Ä–æ–∏–º —Å–ø–ª–∞–π–Ω, —Å–æ–æ—Ç–≤–µ—Ç—Å—Ç–≤—É—é—â–∏–π —Ç–∞–±–ª–∏—á–Ω–æ –∑–∞–¥–∞–Ω–Ω–æ–º—É —Ç–æ–∫—É
	xa1.setcontent(size_kink, xmas1);
	ya1.setcontent(size_kink, ymas1);
	spline1dbuildcubic(xa1, ya1, s1);
	cout << "Done interpolation" << endl;
	cout << "Begin solve wave equation" << endl;
	//–ù–∞—á–∞–ª—å–Ω—ã–µ —É—Å–ª–æ–≤–∏—è 
	hx = (xmax-xmin)/(nx-1);
	ht = hx/2;
	nt = (int)floor((tmax-tmin)/ht)+1; // –û–±—â–µ–µ —á–∏—Å–ª–æ —É–∑–ª–æ–≤ –ø–æ –≤—Ä–µ–º–µ–Ω–∏
	f = (double**)calloc(nx, sizeof(double));//–ß–∏—Å–ª–æ —ç–ª–µ–º–µ–Ω—Ç–æ–≤ –ø–æ –ø–µ—Ä–≤–æ–º—É –∏–Ω–¥–µ–∫—Å—É
	for(int i = 0; i < nx; i++)
	{
		f[i] = (double*)calloc(masnt, sizeof(double));//–ß–∏—Å–ª–æ —ç–ª–µ–º–µ–Ω—Ç–æ–≤ –ø–æ –≤—Ç–æ—Ä–æ–º—É –∏–Ω–¥–µ–∫—Å—É
	};
	// –ó–∞–ø–æ–ª–Ω—è–µ–º –Ω–∞—á–∞–ª—å–Ω—ã–µ —É—Å–ª–æ–≤–∏—è
	// –ï—Å–ª–∏ –æ–¥–∏–Ω–æ—á–Ω—ã–π –∫–∏–Ω–∫
	if (v2 == 0.0)
		for (int x = 0; x < nx; x++)
		{
			f[x][0] = f0_single_kink(xmin+x*hx, tmin, v1, x10, b);
			f[x][1] = f0_single_kink(xmin+x*hx, tmin+ht, v1, x10, b);
		}
	else
		// –ï—Å–ª–∏ –∫–∏–Ω–∫ –¥–≤–æ–π–Ω–æ–π
		for (int x = 0; x < nx; x++)
		{
			f[x][0] = f0_double_kink(xmin+x*hx, tmin, v1, v2, x10, x20, b);
			f[x][1] = f0_double_kink(xmin+x*hx, tmin+ht, v1, v2, x10, x20, b);
		};
	omp_set_num_threads(threads);
	// –í–∫–ª—é—á–∞–µ–º –≤—Ä–µ–º—è
	#pragma omp parallel private (myid) shared (f, threads, masnt, nx, divt, divx, nt, mu, xb, l, alpha_hf)
	{
		double mul = 0.0;
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
					mul = (1 - exp(-(tmin + (t-1)*ht)/tau));
					f[x][t] = 1.0/(1+betta*ht/2)*((ht*ht)*(f[x-1][t-1]+f[x+1][t-1]-2*f[x][t-1])/(hx*hx)+2*f[x][t-1]-f[x][t-2]
<<<<<<< HEAD
							+ betta*ht*f[x][t-2]/2
							- (1 + delta_barrier(xmin + x*hx))*ht*ht*spline1dcalc(s1, f[x][t-1]));//spline1dcalc - —Ç–æ–∫
=======
							+ betta*ht*f[x][t-2]/2 + 0*mul*ht*ht*a*w*w*sin(w*(tmin + (t-1)*ht)) - mul*ht*ht*betta*a*w*cos(w*(tmin + (t-1)*ht))
							- (1 + delta_barrier(xmin + x*hx))*ht*ht*current(f[x][t-1] + mul*a*sin(w*(tmin + (t-1)*ht)))); //spline1dcalc(s1, f[x][t-1]));//spline1dcalc - ÚÓÍ
>>>>>>> 3ab42d3d87c37acc9d9821708d293f81c0f301c4
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
					out_file.open(argv[2]);
				else
					out_file.open(argv[2], std::ios_base::app);
				for(int t = 1; t < masnt; t+=divt)
				{
					for(int x = 0; x < nx - 1; x+=divx)
					{
						out_file << (f[x][t] - f[x][t - 1])/ht << " ";
					};
					out_file << (f[nx - 1][t] - f[nx - 1][t - 1])/ht << endl;
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
	//cin.get();
	return 0;
}

