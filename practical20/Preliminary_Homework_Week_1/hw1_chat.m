% Ex 1
%close all
clear all
clc

% initial parameters definition
t = 0:0.001:10; % time

x0 = 15e-3;                % amplitude of linear translation [m]
a  = 1.5e-3;               % radius of the sphere            [m]

p0 = 1.225;                % density of air                  [kg/m^3] - ToDo add reference
pp = 31;                   % density of styrene              [kg/m^3]
c0 = 346;                  % speed of sound in air           [m/s]
cp = (2350+1120+1840)/3;   % speed of sound in styrene       [m/s]  - https://www.engineeringtoolbox.com/sound-speed-solids-d_713.html
mu0 = 1.81e-5;             % dynamic viscosity of air        [kg/(m·s)] - ToDo add reference
m = (4/3) * pi * a^3 * pp; % mass of the sphere              [kg]

% f1 and f2 computation
f1 = 1 - (p0 * c0^2)/(p0*cp^2);

f2 = 2 * ((pp - p0) / (2 * pp + p0));


% drag force and inertial load computation for 2.5Hz
f=2.5;
x = x0 * sin(2*pi*f*t);

x_dot = x0 * 2 * pi * f * cos(2*pi*f*t);
x_ddot = -x0 * (2*pi*f)^2 * sin(2*pi*f*t);

relative_velocity = x_dot; % Assuming no fluid velocity

Re = abs(relative_velocity) * 2 * a / mu0;

% Using a simplified drag coefficient (assuming laminar flow) for low Re
c_d = 24 ./ Re;

drag = -0.5 .* p0 .* pi .* a.^2 .* abs(relative_velocity).^2 .* c_d;
inertial_load = m * x_ddot;

% plotting
figure(1)
% Plot the first subplot
subplot(2, 1, 1);
plot(x, abs(drag));
xlabel('x [m]');
ylabel('F_d [N]');
title('Drag force along trajectory (f=2.5Hz)');

% Plot the second subplot
subplot(2, 1, 2);
plot(x, abs(inertial_load));
xlabel('x [m]');
ylabel('F_i [N]');
title('Inertial load along trajectory (f=2.5Hz)');
