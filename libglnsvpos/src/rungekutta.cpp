#include <include\libglnsvpos\rungekutta.h>
#include <cmath>
#include <include\libglnsvpos\diffs.h>
#include <iostream>
#include <fstream>

using namespace std;

void rungekutta(double t0, double tn, double *koord0, double *koord){
    double dt = tn-t0;
    double *proizv = new double [6];
    double *K1 = new double [6];
    double *K2 = new double [6];
    double *K3 = new double [6];
    double *K4 = new double [6];
    double *prom1 = new double [6];
    diffs(t0,koord0,proizv);
    for (int j = 0; j<=5; j++){
        K1[j] = dt*proizv[j];
        prom1[j] = koord0[j]+K1[j]/2.0;
    }
    diffs(t0+dt/2.0,prom1,proizv);
    for (int j = 0; j<=5; j++){
        K2[j] = dt*proizv[j];
        prom1[j] = koord0[j]+K2[j]/2.0;
    }
    diffs(t0+dt/2.0,prom1,proizv);
    for (int j = 0; j<=5; j++){
        K3[j] = dt*proizv[j];
        prom1[j] = koord0[j]+K3[j];
    }
    diffs(tn,prom1,proizv);
    for (int j = 0; j<=5; j++){
        K4[j] = dt*proizv[j];
        koord[j] = koord0[j] + (K1[j]+2*K2[j]+2*K3[j]+K4[j])/6.0;
    }
    delete[] K1;
    K1 = nullptr;
    delete[] K2;
    K2 = nullptr;
    delete[] K3;
    K3 = nullptr;
    delete[] K4;
    K4 = nullptr;
    delete[] proizv;
    proizv = nullptr;
    delete[] prom1;
    prom1 = nullptr;
}
