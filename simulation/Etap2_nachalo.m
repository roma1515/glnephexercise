clear all, close all

%Данные из эфемерид спутник номер 14
X = 12234925.29;
Y = -15420709.96;
Z = 16234023.93;

VXt = -683.95519;
VYt = 2115.44323;
VZt = 2525.34008;

AX = -0.0000019;
AY = -0.0000009;
AZ = 0.0000028;

tAY = -4798.9;%ns
gamma = 0;

%Дано
Toe = 13*60*60+45*60+16;% 13 45 16 сек из эфемерид
ti = 12*60*60;%текущее время просмотра
theta_Go =(16*60*60+45*60+16);%из алгоритма расчёта the sidereal time in Greenwich at midnight GMT of a date at which the epoch  is specified
omega_e = 0.7292115*10^-4;% скорость вращения Земли
theta_Ge = theta_Go+omega_e*(Toe - 3*3600);%is the sidereal time at epoch , to which are referred the initial conditions, in Greenwich meridian,10800 сек это три часа

%Пункт 1 переведём в инерциальную систему
Xate = X*cos(theta_Ge) - Y*sin(theta_Ge);
Yate = X*sin(theta_Ge) + Y*cos(theta_Ge);
Zate = Z;
VXate=VXt*cos(theta_Ge)-VYt*sin(theta_Ge)-omega_e*Yate;
VYate=VXt*sin(theta_Ge)+VYt*cos(theta_Ge)+omega_e*Xate;
VZate=VZt;
AX=AX*cos(theta_Ge)-AY*sin(theta_Ge);
AY=AX*sin(theta_Ge)+AY*cos(theta_Ge);
AZ=AZ;
%%%%%%%%%%%%%%%%
%%%%%%
%%%
%Пункт 2 Численное интегрирование дифференциальных уравнений, описывающих движение спутников.
a_e = 6378.136; % Экваториальный радиус Земли (PZ-90)
Mu = 398600.44; %km^3/s^2 гравитационная постоянная (PZ-90)
C20 = -1082.63*10^-6; %Second zonal coefficient of spherical harmonic expression
 
%Расчётные значения
r = sqrt(Xate^2 + Yate^2 + Zate^2);
Mu_strih = Mu/(r^2);
Xate_strih = Xate/r;
Yate_strih = Yate/r;
Zate_strih = Zate/r;
rO = a_e/r;
 
%Дифференциальные уравнения, описывающие движение спутников.
 
dxa_dt=VXate;
dya_dt=VYate;
dza_dt=VZate;
dVXate_dt=-Mu_strih*Xate_strih+1.5*C20*Mu_strih*Xate_strih*(rO^2)*(1-5*Zate_strih^2)+AX;
dVYate_dt=-Mu_strih*Yate_strih+1.5*C20*Mu_strih*Yate_strih*(rO^2)*(1-5*Zate_strih^2)+AY;
dVZate_dt=-Mu_strih*Zate_strih+1.5*C20*Mu_strih*Zate_strih*(rO^2)*(3-5*Zate_strih^2)+AZ;
 
%Метод Рунге-Кутты 
% задать время ->