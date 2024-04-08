clear all
%close all
clc


t = 0.001:0.001:10;
a = 1.5e-3;
x0 = 15e-3;
rho_0 = 1.225;
mu0 = 1.81e-5;             % dynamic viscosity of air        [kg/(m·s)] - ToDo add reference
rho_styrene = 31;
m = 4/3 * pi * a^3 * rho_styrene;

f = 2.5;
x = x0*sin(2*pi*f*t);

x_dot = 2*pi*f*x0*cos(2*pi*f*t);

x_dot_dot = 4*pi*pi*f*f * x0 * (-1) * sin(2*pi*f*t);

Re = abs(x_dot) * rho_0 * 2 * a / mu0;
c_d_1 = (24./Re);

c_d_2 = (((2.6/5)*Re) ./ (1 + (Re./5).^1.52));

c_d_3 = ((0.411*(Re./2.63e5).^-7.94) ./ (1 + (Re./2.63e5).^-8 ));

c_d_4 = ((0.25e-6*Re) ./ (1 + Re./1e6));

c_d = c_d_1 + c_d_2 + c_d_3 + c_d_4;

Fd = - 0.5 * rho_0 * pi * a^2 * c_d .* abs(x_dot).^2 .* Re;

Fin = m * x_dot_dot;

figure(3)
plot(x,abs(Fd))
figure(4)
plot(x,abs(Fin))