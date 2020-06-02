#include <cmath>
#include <include\libglnsvpos\diffs.h>

void diffs(double t, double *koord, double *proizv)
{
    double mu = 3.986004418E+14; // конствнтва гравитационного поля Земли
    double Rz = 6378136; // экваториальный радиус Земли
    double we = 7.292115E-05; // угловая скорость вращения Земли
    double C20 = -1082.62575E-06;
    double t_G0 = (9*3600+18*60+10.5009+3*3600);
    double t_G = t_G0 + we*(t- 3*3600);
    // Ускорения, принимаем постоянными на всем интервале расчета
    double Ax = .186264514923E-05;
    double Ay = .931322574615E-06;
    double Az = -.279396772385E-05;

    double Jsum_x = Ax*cos(t_G)-Ay*sin(t_G);
    double Jsum_y = Ax*sin(t_G)+Ay*cos(t_G);
    double Jsum_z = Az;

    proizv[0] = koord[3];
    proizv[1] = koord[4];
    proizv[2] = koord[5];
    double r = sqrt(koord[0]*koord[0]+koord[1]*koord[1]+koord[2]*koord[2]);
    double mu_ = mu/(r*r);
    double x_ = koord[0]/r;
    double y_ = koord[1]/r;
    double z_ = koord[2]/r;
    double ro = Rz/r;
    proizv[3] = -mu_*x_ + 3/2*C20*mu_*x_*ro*ro*(1-5*z_*z_) + Jsum_x;
    proizv[4] = -mu_*y_ + 3/2*C20*mu_*y_*ro*ro*(1-5*z_*z_) + Jsum_y;
    proizv[5] = -mu_*z_ + 3/2*C20*mu_*z_*ro*ro*(3-5*z_*z_) + Jsum_z;
}
