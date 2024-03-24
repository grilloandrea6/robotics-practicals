% Ex 1
close all
clear all
t = 0:0.01:2*pi;
x0 = 0.015; % m

f = 1;
a = 0.0015; % m

p0 = 1.293; % kg/m^3
pp = 28;    % kg/m^3
mu0 = 1.81e-5; % kg/(m·s) dynamic viscosity of air

x = x0 * sin(2*pi*f*t);
y = 0;



x_dot = x0 * 2 * pi * f * cos(2*pi*f*t);

q_dot_norm = abs(x_dot);

Re = q_dot_norm * p0 * 2 * a / mu0;

c0 = 346;
cp = 1234; %??
cd = 24./Re + 2.6*(Re./5) ./ (1 + (Re./5).^1.52) + 0.411*(Re./2.63e5).^-7.94 ./ (1 + (Re./2.63e5).^-8 ) + 0.25*(Re./1e6) ./ (1 + Re./1e6);


drag = -0.5 .* p0 .* pi .* a.^2 .* q_dot_norm.^2 .* cd .* Re


f1 = 1 - (p0 * c0^2)/(p0*cp^2)

f2 = 2 * ((pp - p0) / (2 * pp + p0))

U = 2 .* pi .* a.^3 .* ()

figure(1)
plot(t, drag)
figure(2)
plot(x, drag)