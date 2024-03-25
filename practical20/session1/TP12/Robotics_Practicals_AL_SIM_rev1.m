%% Robotics Practicals - Numerical simulation
% Controlling the position of the levitated particle by controling the 
% position and oriantation of the levitatotr.
%
% This script enables to control the harmonical motion of the rotational
% motion (motor 1) and translational motion (motor 2). That is, the
% amplitude, frequency and phase for each motor can be tuned.
%
% There are few assumptions regarding the parameters related to the system
% that can be adjusted, such as the particle's density, radius, the gravity
% coefficient and the applied voltage. (if you cahnge the gravity and
% density, make sure you change it throughout the sscript).
%
% After the initialization of the script, first, a static simulation is 
% done to esrimate the intial conditions for the dynamic simulation. In
% that section "Simulate Static equillibrium"  you can set the initial
% angle and position of the levitator. During the simulation Numerical
% Damping is added so that the partticle reach an equilibrium position
% faster. (Set the amplitude and phase for motor 1)
%
% In the following section "Simulate dynamic case - start from the 
% equilibrium point", the dynamic case is simulated. (Set the frequency
% for motor 1, set the amplitude, phase and frequency for motor 2)
%
% After both simulations, figures contatining data will display the results
% including an animation of the simulated case if ANIMATE = 1;.

% Initiate the script
clc; clear; close all

% Define graphics parameters
Fsize = 12;                     % Font Size
FPos = [0.1 0.1 0.8 0.8];       % Figure position
LW = 2;                         % Line Width
ANIMATE = 1;                    % Animate the dynamic response
DS = 1;                         % DownSample for animation visualization (Higher INTEGER numbers lead to faster animations)
% Define parameteres 
V = 5;                          % V - Transducers' voltage amplitude (0-10)
a = 1.5e-3;                     % m - particle radius
rho_sty = 31;                   % The density range is about 28–34 kg/m3: https://en.wikipedia.org/wiki/Polystyrene#:~:text=The%20density%20range%20is%20about%2028%E2%80%9334%20kg%2Fm3.
g = 9.8;                        % m s-2 - Gravity constant
mass = 4/3*pi*a^3*rho_sty;      % Particle mass
Mg = mass*g;                    % Particle weight

% Load simulated forces - used for the simulations
load("Al_Gorkov_1m_1V.mat");    % Assuming 1V, radius = 1m 'ZZgrk','RRgrk','Frgrk','Fzgrk'

% Request prescalar value
% prompt = {'Enter the Prescaler value'};
% dlgtitle = 'Prescaler - to remove the dialog box comment lines 49-52 & uncoment line 53';
% dims = [1 100]; definput = {'1'}; answer = inputdlg(prompt,dlgtitle,dims,definput);
% Prescaler = str2double(answer{1});
Prescaler = 1.97163;
Frgrk = Prescaler*Frgrk;
Fzgrk = Prescaler*Fzgrk;

%% Simulate Static equillibrium
% Set the harmonic signals parameteres - Do not change lines 60-64.
% Translation u =  Amp_u*sin((fex_u*2*pi)*tt+phi_u); (motor 2)
fex_u = 0;                  % Hz - translational oscillation frequency
Amp_u = 0;                  % m - translational oscillation amplitude
phi_u = 0;                  % rad - translational phase
% Rotation th =  Amp_th*sin((fex_th*2*pi)*tt+phi_th); (motor 2)
fex_th = 0;                 % Hz - rotational oscillation frequency
Amp_th = 45/180*pi;       % rad - rotational oscillation amplitude
phi_th = 0;        % rad - rotational phase

w_u = 2*pi*fex_u;           % Hz -> rad/s
w_th = 2*pi*fex_th;         % Hz -> rad/s
% Initial conditions 
IC = [0 0 0 0]';            % [r z dr dz] - initial conditions (ode45)
tt = linspace(0,10,1e3)';    % s - define the time vector (ode45)
opts = odeset('RelTol',1e-6,'AbsTol',1e-8); 

% Simulate 
cN = 1e-4;%3*mass*2*1e-2*(50*2*pi);        % Add Numerical damping
[~,x_0] = ode45(@(t,x) Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk,Fzgrk,cN), tt, IC, opts);

r = x_0(:,1);                       % m - get the r coordinate
z = x_0(:,2);                       % m - get the z coordinate
u = Amp_u*sin(w_u*tt+phi_u);        % m - get the u coordinate
th = Amp_th*sin(w_th*tt+phi_th);    % rad - get the theta coordinate
X = u + r.*cos(th)-z.*sin(th);      % m - compute the particle x position in th
Y = r.*sin(th)+z.*cos(th);          % m - compute the particle y position in th
% Get the particle equilibrium position and use it as the intial condition
% in the dynamic simulation
r0 = x_0(end,1);                    
z0 = x_0(end,2);

% Plot the results of the static simulation
figure(1); clf
subplot(2,2,1)
    plot(tt,X*1e3,'LineWidth',LW);
    ylabel('$x$ (mm)','Interpreter','latex')
    set(gca,'FontSize',Fsize)
subplot(2,2,2)
    plot(tt,Y*1e3,'LineWidth',LW);
    ylabel('$y$ (mm)','Interpreter','latex')
    set(gca,'FontSize',Fsize)
subplot(2,2,3)
    plot(tt,r*1e3,'LineWidth',LW);
    ylabel('$r$ (mm)','Interpreter','latex')
    xlabel('$t$ (s)','Interpreter','latex')
    set(gca,'FontSize',Fsize)
subplot(2,2,4)
    plot(tt,z*1e3,'LineWidth',LW);
    ylabel('$z$ (mm)','Interpreter','latex')
    xlabel('$t$ (s)','Interpreter','latex')
    set(gca,'FontSize',Fsize)
set(gcf,'color','w','units','normalized','position',FPos)
%% Simulate dynamic case - start from the equilibrium point

% % my_ampl = [5e-3 10e-3 15e-3];
% my_freq = [1 2 3 4 5];
my_ampl = 15e-3;
my_freq = 5;
[my_ampl, my_freq] = ndgrid(my_ampl, my_freq);

for my_idx = 1:length(my_ampl(:))
    % Translation u =  Amp_u*sin((fex_u*2*pi)*tt+phi_u); (motor 2)
%     fex_u = 5;                  % Hz - translational oscillation frequency
%     Amp_u = 10e-3;              % m - translational oscillation amplitude
%     phi_u = -pi/2;         % rad - translational phase
%     % Rotation th =  Amp_th*sin((fex_th*2*pi)*tt+phi_th); (motor 1)
%     fex_th = 0;     % Hz - rotational oscillation frequency

    % Translation u =  Amp_u*sin((fex_u*2*pi)*tt+phi_u); (motor 2)
    fex_u = my_freq(my_idx);                  % Hz - translational oscillation frequency
    Amp_u = my_ampl(my_idx);              % m - translational oscillation amplitude
    phi_u = -pi/2;         % rad - translational phase
    % Rotation th =  Amp_th*sin((fex_th*2*pi)*tt+phi_th); (motor 1)
    Amp_th = 60/180*pi;
    fex_th = fex_u;                 % Hz - rotational oscillation frequency
    phi_th = -pi/2;


    w_u = 2*pi*fex_u;           % Hz -> rad/s
    w_th = 2*pi*fex_th;         % Hz -> rad/s
    % Automatically set the time vector according to the frequencies
    Maxf = max([fex_u;fex_th]); 
    minf = min([fex_u;fex_th]);
    if ~minf
        minf = Maxf;
    end
    % Initial conditions 
    IC = [r0 z0 0 0]';          % [r z dr dz] - initial conditions (ode45)
    dt = 1e-2/Maxf;
    tt = (0:dt:30/minf)';       % s - define the time vector (ode45)
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    
    % Simulate 
    cN = 0;						% remove Numerical damping
    [~,x_0] = ode45(@(t,x) Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk,Fzgrk,cN), tt, IC, opts);
    
    r = x_0(:,1);                       % m - get the r coordinate
    z = x_0(:,2);                       % m - get the z coordinate
    u = Amp_u*sin(w_u*tt+phi_u);        % m - get the u coordinate
    th = Amp_th*sin(w_th*tt+phi_th);    % rad - get the theta coordinate
    X = u + r.*cos(th)-z.*sin(th);      % m - compute the particle x position in th
    Y = r.*sin(th)+z.*cos(th);          % m - compute the particle y position in th
    
    % Plot theresults of the dynamic symulation
    figure(2); clf
    subplot(2,2,1)
        plot(tt,X*1e3,'LineWidth',LW);
        ylabel('$x$ (mm)','Interpreter','latex')
        set(gca,'FontSize',Fsize)
        xlim([tt(1) tt(end)])
    subplot(2,2,2)
        plot(tt,Y*1e3,'LineWidth',LW);
        ylabel('$y$ (mm)','Interpreter','latex')
        set(gca,'FontSize',Fsize)
        xlim([tt(1) tt(end)])
    subplot(2,2,3)
        plot(tt,r*1e3,'LineWidth',LW);
        ylabel('$r$ (mm)','Interpreter','latex')
        xlabel('$t$ (s)','Interpreter','latex')
        set(gca,'FontSize',Fsize)
        xlim([tt(1) tt(end)])
    subplot(2,2,4)
        plot(tt,z*1e3,'LineWidth',LW);
        ylabel('$z$ (mm)','Interpreter','latex')
        xlabel('$t$ (s)','Interpreter','latex')
        set(gca,'FontSize',Fsize)
        xlim([tt(1) tt(end)])
    set(gcf,'color','w','units','normalized','position',FPos)

     filename = "plot5V-5Hz-15mm-freqth5Hz-ampth60-ph-90"; %sprintf('plotV%d-f%d-u%d-fth%d-ampth%d-phth%d', round(V*1000),round(my_freq(my_idx)*10^4),round(1000*my_ampl(my_idx),round(1000*fex_th,roundTies="toEven"),round(1000*Amp_th,roundTies="toEven"),round(1000*phi_th,roundTies="toEven")));
     saveas(gcf,filename,'eps')
     saveas(gcf,filename,'png')
     


    
    % Plot the steady state - last slow cycles or common cycle
    cycles = 1;
    Tss = lcm(round(1e3/minf),round(1e3/Maxf));
    Tss = Tss*1e-3;
    if Tss>tt(end)
        Tss = cycles/minf;
    end
    tt_ss = tt(tt >= tt(end) - Tss);
    X_ss = X(tt >= tt(end) - Tss);
    Y_ss = Y(tt >= tt(end) - Tss);
    u_ss = u(tt >= tt(end) - Tss);
    
    figure(3); clf
    subplot(2,2,1)
        plot(X_ss*1e3,Y_ss*1e3,'LineWidth',LW);
        xlabel('$x$ (mm)','Interpreter','latex')
        ylabel('$y$ (mm)','Interpreter','latex')
        set(gca,'FontSize',Fsize)
    subplot(2,2,2)
        plot((X_ss-u_ss)*1e3,Y_ss*1e3,'LineWidth',LW);
        ylabel('$y$ (mm)','Interpreter','latex')
        xlabel('$x-u$ (mm)','Interpreter','latex')
        set(gca,'FontSize',Fsize)
    subplot(2,2,3)
        plot(X_ss*1e3,Y_ss*1e3,'LineWidth',LW);
        ylabel('$y$ (mm)','Interpreter','latex')
        xlabel('$x$ (mm)','Interpreter','latex')
        axis equal
        set(gca,'FontSize',Fsize)
    subplot(2,2,4)
        plot((X_ss-u_ss)*1e3,Y_ss*1e3,'LineWidth',LW);
        ylabel('$y$ (mm)','Interpreter','latex')
        xlabel('$x-u$ (mm)','Interpreter','latex')
        axis equal
        set(gca,'FontSize',Fsize)
    set(gcf,'color','w','units','normalized','position',FPos)
end

%% Animate the results
if ANIMATE
    XML = max(abs(u(:)))*1e3 + 10;
    YML = 13;
    LT = 14.8e-3;                           % Distance between the transducers
    DT = 10e-3;                             % The diameter of each transducer
    ind_anim = 1:DS:length(tt);             % Downsample the data
    hf = figure(4); clf
    set(hf,'units','Normalized','position',[0.1 0.1 0.8 0.8],'color','w')
    
    
    for kk = 1:length(ind_anim)
        ind = ind_anim(kk);
        % Plot The particle
        plot(1e3*X(ind),1e3*Y(ind),'o','MarkerEdgeColor','blue','MarkerFaceColor','b','MarkerSize',6); hold all

        % Plot the transducers
        plot(1e3*[u(ind)+0.5*LT*sin(th(ind));u(ind)-0.5*LT*sin(th(ind))],1e3*[-0.5*LT*cos(th(ind));0.5*LT*cos(th(ind))],'--k');

        plot(1e3*[u(ind)-0.25*LT*cos(th(ind));u(ind)+0.25*LT*cos(th(ind))],1e3*[-0.25*LT*sin(th(ind));0.25*LT*sin(th(ind))],'--k');

        plot(1e3*[u(ind)+0.5*LT*sin(th(ind))-0.5*DT*cos(th(ind));u(ind)+0.5*LT*sin(th(ind))+0.5*DT*cos(th(ind))],...
             1e3*[-0.5*LT*cos(th(ind))-0.5*DT*sin(th(ind));-0.5*LT*cos(th(ind))+0.5*DT*sin(th(ind))],'k','LineWidth',3);
        plot(1e3*[u(ind)-0.5*LT*sin(th(ind))-0.5*DT*cos(th(ind));u(ind)-0.5*LT*sin(th(ind))+0.5*DT*cos(th(ind))],...
             1e3*[0.5*LT*cos(th(ind))-0.5*DT*sin(th(ind));0.5*LT*cos(th(ind))+0.5*DT*sin(th(ind))],'k','LineWidth',3);
        % Highlight the center
        plot(1e3*u(ind),0,'o','MarkerEdgeColor','red','MarkerSize',20); 
                 
        xline(0,'--','color',[0.8 0.8 0.8])
        yline(0,'--','color',[0.8 0.8 0.8])

        xlabel('x (mm)')
        ylabel('y (mm)')
        title(sprintf('t = %.3f s',tt(ind)))
        axis equal
        axis([-XML XML -YML YML])    
        drawnow; shg
        hold off
%         pause(0.03)
    end     
end

%% Functions
function dxdt = Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk,Fzgrk,cN)

dxdt = 0*x;
r = x(1);
z = x(2);
dr = x(3);
dz = x(4);

% Constants
rho_sty = 31;       % The density range is about 28–34 kg/m3: https://en.wikipedia.org/wiki/Polystyrene#:~:text=The%20density%20range%20is%20about%2028%E2%80%9334%20kg%2Fm3.
rho_air = 1.201328767349156;    % kg/m3    amb2prop(pa=101325,temp=20,HHr=30,f_ex=40e3)
mu_air = 1.813508876550286e-05; % N s/m2   amb2prop(pa=101325,temp=20,HHr=30,f_ex=40e3)
g = 9.8;
m = 4/3*pi*a^3*rho_sty;

% Compute the inputes
% u = Amp_u*sin(w*t+phi_u);
du =  w_u*Amp_u*cos(w_u*t+phi_u);
ddu =  -w_u^2*Amp_u*sin(w_u*t+phi_u);

th = Amp_th*sin(w_th*t+phi_th);
dth =  w_th*Amp_th*cos(w_th*t+phi_th);
ddth =  -w_th^2*Amp_th*sin(w_th*t+phi_th);

% Gorkov's contribution
FrN = interp2(RRgrk,ZZgrk,Frgrk,r,z);
FzN = interp2(RRgrk,ZZgrk,Fzgrk,r,z);

FrN = a^3*V^2*FrN;
FzN = a^3*V^2*FzN;

% Compute the velocity and Re
v = (((-1).*du+((-1).*dr+dth.*z).*cos(th)+(dz+dth.*r).*sin(th)).^2+(( ...
  dz+dth.*r).*cos(th)+(dr+(-1).*dth.*z).*sin(th)).^2).^(1/2);
Re = v*rho_air*(2*a)/mu_air;

% Compute Cd: https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdfi
if Re~=0
    cd = 0.444929E43.*(1+0.2289E44.*Re.^(-8)).^(-1).*Re.^(-0.794E1)+24.* ...
      Re.^(-1)+0.25E-6.*(1+(1/1000000).*Re).^(-1).*Re+0.52E0.*Re.*(1+ ...
      0.866095E-1.*Re.^0.152E1).^(-1);
else 
    cd = 0;
end

% computes the derivative of the states z
dxdt(1) = x(3);     % dr/dt
dxdt(2) = x(4);     % dz/dt

dxdt(3) = (-1/2).*m.^(-1).*(a.^2.*cd.*pi.*rho_air.*v.*(dr+(-1).*dth.*z+du.* ...
  cos(th))+2.*((-1).*FrN+m.*((-1).*dth.*(2.*dz+dth.*r)+(-1).* ...
  ddth.*z+ddu.*cos(th)+g.*sin(th)))) - cN/m*dr;
dxdt(4) = (-1/2).*m.^(-1).*(a.^2.*cd.*pi.*rho_air.*v.*(dz+dth.*r+(-1).*du.* ...
  sin(th))+2.*((-1).*FzN+m.*(2.*dr.*dth+ddth.*r+(-1).*dth.^2.*z+ ...
  g.*cos(th)+(-1).*ddu.*sin(th)))) - cN/m*dz;
end

