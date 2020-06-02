#include <include\libglnsvpos\glnsvpos.h>
#include <include\libglnsvpos\rungekutta.h>

#include <iostream>
#include <cmath>
#include <ostream>

using namespace std;

void glns_coord(double **koord_n)
{
    // 14 20  2 10 13 45  0.0  .479789450765E-04  .000000000000E+00  .495000000000E+05
    //      .122349252930E+05 -.683955192566E+00  .186264514923E-08  .000000000000E+00
    //     -.154207099609E+05  .211544322968E+01  .931322574615E-09 -.700000000000E+01
    //      .162340239258E+05  .252534008026E+01 -.279396772385E-08  .000000000000E+00
    double t0 = 49500 + 18 + 3*3600;
    double ts = 12*3600 + 3*3600;
    double te = ts + 12*3600;
    double dt = 1/10.0; //шаг по времени
    double we = 7.292115E-05; // угловая скорость вращения Земли
    int n = (int) abs((te-ts)/dt);
    double *t = new double [n];
    double *koord_p = new double [6];
    for (int i = 0; i < n; i++)
    {
        t[i] = ts + i*dt;
    }
    int n_eph = (int) (t0-ts)/dt;
    // начальные условия (коорд и скорости)
    koord_p[0] = 0.122349252930E+08;
    koord_p[1] = -0.154207099609E+08;
    koord_p[2] = 0.162340239258E+08;
    koord_p[3] = -0.683955192566E+03;
    koord_p[4] = 0.211544322968E+04;
    koord_p[5] = 0.252534008026E+04;

    double t_G0 = (9*3600+18*60+10.5009+3*3600); //9:18:10.5009; Истинное звездное время на гринвичевскую полночь текущей даты
    double t_G = t_G0 + we*(t[n_eph-1]- 3*3600);
    double cosTg = cos(t_G);
    double sinTg = sin(t_G);
    // пересчет координат из ПЗ-90 в ECI
    koord_n[n_eph-1][0] = koord_p[0]*cosTg - koord_p[1]*sinTg;
    koord_n[n_eph-1][1] = koord_p[0]*sinTg + koord_p[1]*cosTg;
    koord_n[n_eph-1][2] = koord_p[2];
    koord_n[n_eph-1][3] = koord_p[3]*cosTg - koord_p[4]*sinTg - we*koord_n[n_eph-1][1];
    koord_n[n_eph-1][4] = koord_p[3]*sinTg + koord_p[4]*cosTg + we*koord_n[n_eph-1][0];
    koord_n[n_eph-1][5] = koord_p[5];
    delete[] koord_p;
    koord_p = nullptr;
    //rungekutta(double t0, double tn, double *koord0, double *koord)
    for (int i = n_eph-1; i > 0; i--)
    {
        rungekutta(t[i], t[i-1], koord_n[i], koord_n[i-1]);
    }
    for (int i = n_eph-1; i < n-1 ; i++)
    {
        rungekutta(t[i], t[i+1], koord_n[i], koord_n[i+1]);
    }
    delete[] t;
    t = nullptr;
    delete[] koord_p;
    koord_p = nullptr;
}
