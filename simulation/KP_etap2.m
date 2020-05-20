% 14 20  2 10 13 45  0.0  .479789450765E-04  .000000000000E+00  .495000000000E+05
%      .122349252930E+05 -.683955192566E+00  .186264514923E-08  .000000000000E+00
%     -.154207099609E+05  .211544322968E+01  .931322574615E-09 -.700000000000E+01
%      .162340239258E+05  .252534008026E+01 -.279396772385E-08  .000000000000E+00

clear all;
close all;

tic;
%������������ ���������
del_t = 1; % ��� ��������
mu = 3.986004418E+14; % ��������� ��������������� ���� �����
we = 7.292115E-05; % ������� �������� �������� �����
Rz = 6371000; % ������ �����
C20 = -1082.62575E-06; % ����������� ��� ������ ��������� ��������� ����������
%������������� � ��� �� ����������� ��������

%���������
NS = 14; % ����� ��������
T_Omega = .495000000000E+05 + 18 + 3*3600; % ����� ������� ��������
X = .122349252930E+08; 
Y = -.154207099609E+08;
Z = .162340239258E+08;
Vx = -.683955192566E+03;
Vy = .211544322968E+04;
Vz = .252534008026E+04;
Ax = .186264514923E-05;
Ay = .931322574615E-06;
Az = -.279396772385E-05;


num = fix(12*3600/del_t); % ���������� �������� �� 12 ������� �������� �������
num_eph = fix((T_Omega-12*3600-3*3600)/del_t); % ����� ������� � ������������ �������
tt = del_t.*(1:1:num); %������ �������� ������� 
tt = tt + 12*60*60+3*60*60; %������ �������� ������� 
koord_PZ = zeros(num,3); 
koord_ECI = zeros(num,3);
LL_potr = [55.756735, 37.703177 170]; % ���������� �����������

t_G0 = (9*3600+18*60+10.5009+3*3600); %9:18:10.5009; �������� �������� ����� �� ������������� ������� ������� ����
t_G = t_G0 + we*(tt(num_eph)- 3*3600);

% �������� ��������� �� ��-90 � ECI
Xa = X*cos(t_G) - Y*sin(t_G); 
Ya = X*sin(t_G) + Y*cos(t_G);
Za = Z;
Vxa = Vx*cos(t_G) - Vy*sin(t_G) - we*Ya;
Vya = Vx*sin(t_G) + Vy*cos(t_G) + we*Xa;
Vza = Vz;
Yn = [Xa Ya Za Vxa Vya Vza];

% �������������� ������� �����-�����
[t, Yn1] = ode45('proizv', T_Omega:-del_t:tt(1), Yn);
koord_ECI(1:num_eph,:) = Yn1(end:-1:1,1:3);
[t, Yn1] = ode45('proizv', T_Omega:del_t:tt(num), Yn);
koord_ECI(num_eph:end,:) = Yn1(1:end,1:3);

% �������� ���������� ��������� �� ECI � ��-90
for i = 1:num
    t_G = t_G0 + we*(tt(i)- 3*3600);
    koord_PZ(i,1) = koord_ECI(i,1)*cos(t_G) + koord_ECI(i,2)*sin(t_G);
    koord_PZ(i,2) = -koord_ECI(i,1)*sin(t_G) + koord_ECI(i,2)*cos(t_G);
    koord_PZ(i,3) = koord_ECI(i,3);
end

% �������� ��������� �� ��-90 � WGS84 (��� ��������� SkyView)
ppb = 1e-9;
mas = 1e-3/206264.8; % [���]
M_WGS84 = [-3*ppb -353*mas -4*mas;
    353*mas -3*ppb 19*mas;
    4*mas -19*mas -3*ppb];

koord_WGS84 = koord_PZ.'; % ������� � �������-�������
for i = 1:length(koord_WGS84(1,:))
    koord_WGS84(:,i) =  koord_WGS84(:,i) + M_WGS84 * koord_WGS84(:,i) + [0.07; -0; -0.77];
end
koord_WGS84 = koord_WGS84.'; % ������� � �������-������

% �������� ��������� �� WGS84 � SkyView
X = zeros(num,1);
Y = zeros(num,1);
Z = zeros(num,1);
r = zeros(num,1);
teta = zeros(num,1);
phi = zeros(num,1);
for i = 1:length(koord_WGS84(:,1))
    [X(i), Y(i), Z(i)] = ecef2enu(koord_WGS84(i,1),koord_WGS84(i,2),koord_WGS84(i,3),LL_potr(1),LL_potr(2),LL_potr(3),wgs84Ellipsoid);
    if Z(i) > 0
        r(i) = sqrt(X(i)^2 + Y(i)^2 + Z(i)^2);
        teta(i) = acos(Z(i)/r(i));
        if X(i) > 0
            phi(i) = -atan(Y(i)/X(i))+pi/2;
        elseif (X(i)<0)&&(Y(i)>0)
            phi(i) = -atan(Y(i)/X(i))+3*pi/2;
        elseif (X(i)<0)&&(Y(i)<0)
            phi(i) = -atan(Y(i)/X(i))-pi/2;
        end
    else
        teta(i) = NaN;
        r(i) = NaN;
        phi(i) = NaN;
    end
end

[Xsf, Ysf, Zsf] = sphere(25);
alfa1 = pi/180.*(1:359)';
beta1 = 85.*ones(359,1);

plot3(koord_PZ(:,1), koord_PZ(:,2), koord_PZ(:,3));
hold on
grid on
title('��������� �������� � �� ��-90');
xlabel('OX, ��');
ylabel('OY, ��');
zlabel('OZ, ��');
axis('square');
axis('equal');
surf(Xsf.*Rz, Ysf*Rz, Zsf*Rz);
hold off

figure;
plot3(koord_ECI(:,1), koord_ECI(:,2), koord_ECI(:,3));
hold on
grid on
title('��������� �������� � �� ECI');
xlabel('OX, ��');
ylabel('OY, ��');
zlabel('OZ, ��');
axis('square');
axis('equal');
surf(Xsf.*Rz, Ysf*Rz, Zsf*Rz);
hold off

figure;
ax = polaraxes;
hold on
polarplot(ax,phi,teta*180/pi,'r')
polarplot(ax,alfa1,beta1,'b')
hold off
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
title('��������� �������� SkyView ������� �14')

toc;