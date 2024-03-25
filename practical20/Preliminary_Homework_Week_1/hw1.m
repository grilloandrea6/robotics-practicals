% Ex 1
close all
clear all

% initial parameters definition
t = 0:0.01:10; % time

x0 = 0.015;                % amplitude of linear translation [m]
a = 0.0015;                % radius of the sphere            [m]

p0 = 1.293;                % density of air                  [kg/m^3] - ToDo add reference
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


Re = abs(x_dot) * p0 * 2 * a / mu0;

cd = 24./Re + 2.6*(Re./5) ./ (1 + (Re./5).^1.52) + 0.411*(Re./2.63e5).^-7.94 ./ (1 + (Re./2.63e5).^-8 ) + 0.25*(Re./1e6) ./ (1 + Re./1e6);

drag = -0.5 .* p0 .* pi .* a.^2 .* abs(x_dot).^2 .* cd .* Re;
inertial_load = m * x_ddot;

% plotting
figure(1)
% Plot the first subplot
subplot(2, 1, 1);
plot(x, abs(drag));
xlabel('x [mm]');
ylabel('F_d [N]');
title('Drag force along trajectory (f=2.5Hz)');

% Plot the second subplot
subplot(2, 1, 2);
plot(x, abs(inertial_load));
xlabel('x [mm]');
ylabel('F_i [N]');
title('Inertial load along trajectory (f=2.5Hz)');


% maximum drag force and inertial load for range of frequencies

inload_range = [];
drag_range = [];
frequencies = linspace(0.1,5,50); % Hz

for i = 1 : length(frequencies)
    f = frequencies(i);
    
    x = x0 * sin(2*pi*f*t);
    
    x_dot = x0 * 2 * pi * f * cos(2*pi*f*t);
    x_ddot = -x0 * (2*pi*f)^2 * sin(2*pi*f*t);
    
    
    Re = abs(x_dot) * p0 * 2 * a / mu0;
    
    cd = 24./Re + 2.6*(Re./5) ./ (1 + (Re./5).^1.52) + 0.411*(Re./2.63e5).^-7.94 ./ (1 + (Re./2.63e5).^-8 ) + 0.25*(Re./1e6) ./ (1 + Re./1e6);
    
    drag = -0.5 .* p0 .* pi .* a.^2 .* x_dot.^2 .* cd .* Re;
    inertial_load = m * x_ddot;

    inload_range = [inload_range max(abs(inertial_load))];
    drag_range = [drag_range max(abs(drag))];
    
end

figure(2)
subplot(3,1,1)
plot(frequencies,drag_range)
title('Maximum drag force vs frequency')
xlabel('f [Hz]')
ylabel('F_d [N]')

subplot(3,1,2)
plot(frequencies,inload_range)
title('Maximum inertial load vs frequency')
xlabel('f [Hz]')
ylabel('F_i [N]')

subplot(3,1,3)
plot(frequencies,drag_range, "r")
hold on
plot(frequencies,inload_range, "b")

% plot intertial load and drag as a function of frequency
% 3 groups of 2 plots


% intertial leads to escape
% but there is a point around 0.5Hz drag force is higher


% put a reference for each constant we used

