%% Robotics Practicals - Estimate the Gor'kov potenital prescalar  
% The prescalar is estimated using experimentally measured data and 
% numerical simulations. Since in practice the trap does not coincide with 
% the center of rotation, and it is not aligned with the axis of rotation,
% we slightly relax the conditions. To estimate the prescalar we minimize 
% the relative error between the position of the particle as it was 
% estimated experimentaly and comuted numerically (numerical simulation). 
% The relative error is computed relative to the 10V case. 

% For simplicity we define:
% r10S - the simulated r position @ 10V
% z10S - the simulated z position @ 10V
% rS - the simulated r position @ some value of V
% zS - the simulated z position @ some value of V
% r10N - the measured r position @ 10V
% z10N - the measured z position @ 10V
% rN - the measured r position @ some value of V
% zN - the measured z position @ some value of V
% Sim_diff = sqrt( (r10S-rS)^2 + (z10S-zS)^2 ) - The position difference for the simulated case
% Exp_diff = sqrt( (r10N-rN)^2 + (z10N-zN)^2 ) - The position difference for the experimental case
%
% RErr = abs(Exp_diff - Sim_diff)/abs(Sim_diff) - The relative error
% between the simulated and experimental data.

% Initiate the script
clc; clear; close all
COLL = lines(4);
% Load the PNG's processed data file
[File, path] = uigetfile('*.mat','Load processed data files','MultiSelect','Off');
load([path File])

% Load simulated forces - used for the simulations
load("Al_Gorkov_1m_1V.mat");    % Assuming 1V, rad = 1m 'ZZgrk','RRgrk','Frgrk','Fzgrk'

% Prepare estimated data matrix [V r z] and fill it
Data_PNG = [V_V(ind1) r(ind1) z(ind1);
            V_V(ind2) r(ind2) z(ind2);
            V_V(ind3) r(ind3) z(ind3)];
    
%% Use fminbnd to estimate the prescaler that minimize the relative error
% relative particle's position - experiment vs. simulation
[Prescaler_opt,fval,exitflag,output] = fminbnd(@(x) My_opt_fun(x,Data_PNG,RRgrk,ZZgrk,Frgrk,Fzgrk,aN,O_angle_radN,ind1,ind2,ind3),0.1,5);

msg = msgbox({'Optimization Completed'; strcat('Prescaler=',num2str(Prescaler_opt))},'','modal');


%% Plot the results
prescaler = Prescaler_opt;
Frgrk_opt = prescaler*Frgrk;
Fzgrk_opt = prescaler*Fzgrk;

% Simulate Static equillibrium 
fex_u = 0;                  % Hz - oscillation frequency
Amp_u = 0;                  % m - oscillation amplitude
phi_u = 0;                  % rad - phase
% Rotation th =  Amp_th*sin((fex_th*2*pi)*tt+phi_th);
fex_th = 0;                 % Hz - oscillation frequency
phi_th = 90/180*pi;         % rad - phase
w_u = 2*pi*fex_u;
w_th = 2*pi*fex_th;
% Initial conditions 
IC = [0 0 0 0]';          	% [r z dr dz]
cN = 1e-4;	               	% Add Numerical damping
tt = [0 10];              	% s - define the time vector
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
RErr = zeros(length(ind1)+length(ind2)+length(ind3) - (~isempty(ind1) + ~isempty(ind2) + ~isempty(ind3)),1);
rSV = nan + zeros(length(RErr)+3,1);
zSV = rSV;
VsV = rSV;

COL = jet(length(V_V));
figure(3); clf
Sim_ind = 1;
SimT_ind = 1;

for P_ind = 1:(~isempty(ind1) + ~isempty(ind2) + ~isempty(ind3))
    a = aN(P_ind)*1e-3;             % m - the particle radius
    Amp_th = O_angle_radN(P_ind);  	% rad - rotation angle
    
    switch P_ind                    % use the indecies of the particle
        case 1
            ind = ind1;
        case 2
            ind = ind2;
        case 3
            ind = ind3;
    end
    
    % Simulate for 10V and store the result as the reference
    V = Data_PNG(ind(end),1);         % V -tThe voltage used
    
    [t_temp,x_0] = ode45(@(t,x) Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk_opt,Fzgrk_opt,cN), tt, IC, opts);
    figure(2);
    subplot(1,2,1)
        plot(t_temp,x_0(:,1)*1e3);
        xlabel t(s)
        ylabel r(mm)
    subplot(1,2,2)
        plot(t_temp,x_0(:,2)*1e3);
        xlabel t(s)
        ylabel z(mm)
        drawnow; shg
    
    r10S = x_0(end,1)*1e3;      % m-->mm
    z10S = x_0(end,2)*1e3;      % m-->mm
    
    rSV(SimT_ind) = r10S;
    zSV(SimT_ind) = z10S;
    VsV(SimT_ind) = V;
    SimT_ind = SimT_ind+1;

    figure(3);
    subplot(1,2,1)
        plot(VsV,zSV,'--o');
        xlabel V
        ylabel z
    subplot(1,2,2)
        plot(VsV,rSV,'--o');
        xlabel V
        ylabel r
        drawnow; shg

    % Simulate for other voltages
    for kk = 1:length(ind)-1    %[V r z]
        V = Data_PNG(ind(kk),1);
        
        [t_temp,x_0] = ode45(@(t,x) Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk_opt,Fzgrk_opt,cN), tt, IC, opts);
        figure(2);
        subplot(1,2,1)
            plot(t_temp,x_0(:,1)*1e3);
            xlabel t(s)
            ylabel r(mm)
        subplot(1,2,2)
            plot(t_temp,x_0(:,2)*1e3);
            xlabel t(s)
            ylabel z(mm)
             drawnow; shg

        rS = x_0(end,1)*1e3;      % m-->mm
        zS = x_0(end,2)*1e3;      % m-->mm
        
        rSV(SimT_ind) = rS;
        zSV(SimT_ind) = zS;
        VsV(SimT_ind) = V;
        SimT_ind = SimT_ind+1;

        if ~isnan(rS)
        figure(3);
        subplot(1,2,1)
            plot(VsV,zSV,'--o');
            xlabel V
            ylabel z
        subplot(1,2,2)
            plot(VsV,rSV,'--o');
            xlabel V
            ylabel r
            drawnow; shg

        end
        % Compute the differenc and relative error
        Sim_diff = sqrt( (r10S-rS)^2 + (z10S-zS)^2);
        Exp_diff = sqrt( (Data_PNG(ind(end),2)-Data_PNG(ind(kk),2))^2 + (Data_PNG(ind(end),3)-Data_PNG(ind(kk),3))^2);
        RErr(Sim_ind) = abs(Exp_diff - Sim_diff)/abs(Sim_diff);
%                    
%         fprintf('Particle %i, iter %i of %i, Err = %f\n',P_ind,kk,length(ind)-1,RErr(Sim_ind));
        Sim_ind = Sim_ind+1;
    end
end


%%
ZM = z(V_V==10);
RM = r(V_V==10);
FontSize = 10;
figure(100); clf
subplot(1,3,1)
    h1 = plot(V_V(ind1),z(ind1)-ZM(1),V_V(ind2),z(ind2)-ZM(2),V_V(ind3),z(ind3)-ZM(3)); hold all
    h1(4) = plot(VsV([2:10, 1]),zSV([2:10, 1])-zSV(1));
    set(h1(1),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(1,:))
    set(h1(2),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(2,:))
    set(h1(3),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(3,:))
    set(h1(4),'linestyle','--','Marker','p','LineWidth',1,'MarkerFaceColor',COLL(4,:))    
    xlabel('$V$(V)','interpreter','Latex')
    ylabel('$\Delta z$(mm)','interpreter','Latex')
    hl = legend('$a_1$','$a_2$','$a_3$','Sim','Location','Southeast');
    set(hl,'interpreter','Latex','FontSize',FontSize)
    set(gca,'FontSize',FontSize);
subplot(1,3,2)
    h2 = plot(V_V(ind1),r(ind1)-RM(1),V_V(ind2),r(ind2)-RM(2),V_V(ind3),r(ind3)-RM(3),'--o'); hold all
    h2(4) = plot(VsV([2:10, 1]),rSV([2:10, 1])-rSV(1));
    set(h2(1),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(1,:))
    set(h2(2),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(2,:))
    set(h2(3),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(3,:))
    set(h2(4),'linestyle','--','Marker','p','LineWidth',1,'MarkerFaceColor',COLL(4,:))
    xlabel('$V$(V)','interpreter','Latex')
    ylabel('$\Delta r$(mm)','interpreter','Latex')
    set(gca,'FontSize',FontSize);
subplot(1,3,3)
    h3 = plot(V_V(ind1(1:end-1)),RErr(ind1(1:end-1)),...
     V_V(ind2(1:end-1)),RErr(ind2(1:end-1)-1),...
     V_V(ind3(1:end-1)),RErr(ind3(1:end-1)-2));
    yline(mean(RErr),'--r','LineWidth',1)
    set(h3(1),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(1,:))
    set(h3(2),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(2,:))
    set(h3(3),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(3,:))
    xlim([5.5 10])
    xlabel('$V$(V)','interpreter','Latex')
    ylabel('R.E.(\%)','interpreter','Latex')
    set(gca,'FontSize',FontSize);
set(gcf,'color','w','units','centimeters','position',[2 2 20 8])
% exportgraphics(gcf,'Figure7.png','Resolution',300)
figure(200); clf
subplot(1,2,1)
    h1 = plot(V_V(ind1),z(ind1),V_V(ind2),z(ind2),V_V(ind3),z(ind3)); hold all
    h1(4) = plot(VsV([2:10, 1]),zSV([2:10, 1]));
    set(h1(1),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(1,:))
    set(h1(2),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(2,:))
    set(h1(3),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(3,:))
    set(h1(4),'linestyle','--','Marker','p','LineWidth',1,'MarkerFaceColor',COLL(4,:))    
    xlabel('$V$(V)','interpreter','Latex')
    ylabel('$z$(mm)','interpreter','Latex')
    hl = legend('$a_1$','$a_2$','$a_3$','Sim','Location','Southeast');
    set(hl,'interpreter','Latex','FontSize',FontSize)
    set(gca,'FontSize',FontSize);
subplot(1,2,2)
    h2 = plot(V_V(ind1),r(ind1),V_V(ind2),r(ind2),V_V(ind3),r(ind3),'--o'); hold all
    h2(4) = plot(VsV([2:10, 1]),rSV([2:10, 1]));
    set(h2(1),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(1,:))
    set(h2(2),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(2,:))
    set(h2(3),'linestyle','--','Marker','o','LineWidth',1,'MarkerFaceColor',COLL(3,:))
    set(h2(4),'linestyle','--','Marker','p','LineWidth',1,'MarkerFaceColor',COLL(4,:))
    xlabel('$V$(V)','interpreter','Latex')
    ylabel('$r$(mm)','interpreter','Latex')
    set(gca,'FontSize',FontSize);
set(gcf,'color','w','units','centimeters','position',[2 2 20 8])

%% Functions
function Val = My_opt_fun(x,Data_PNG,RRgrk,ZZgrk,Frgrk,Fzgrk,aN,O_angle_radN,ind1,ind2,ind3)
    % This function returns average relative error value for all the cases

    prescaler = x;
    Frgrk = prescaler*Frgrk;
    Fzgrk = prescaler*Fzgrk;
    
    % Simulate Static equillibrium 
    fex_u = 0;                  % Hz - oscillation frequency
    Amp_u = 0;                  % m - oscillation amplitude
    phi_u = 0;                  % rad - phase
    % Rotation th =  Amp_th*sin((fex_th*2*pi)*tt+phi_th);
    fex_th = 0;                 % Hz - oscillation frequency
    phi_th = 90/180*pi;         % rad - phase
    w_u = 2*pi*fex_u;
    w_th = 2*pi*fex_th;
    % Initial conditions 
    IC = [0 0 0 0]';          	% [r z dr dz]
    cN = 1e-4;	               	% Add Numerical damping
    tt = [0 10];              	% s - define the time vector
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    RErr = zeros(length(ind1)+length(ind2)+length(ind3) - (~isempty(ind1) + ~isempty(ind2) + ~isempty(ind3)),1);
   
    Sim_ind = 1;
    for P_ind = 1:(~isempty(ind1) + ~isempty(ind2) + ~isempty(ind3))
        a = aN(P_ind)*1e-3;             % m - the particle radius
        Amp_th = O_angle_radN(P_ind);  	% rad - rotation angle
        
        switch P_ind                    % use the indecies of the particle
            case 1
                ind = ind1;
            case 2
                ind = ind2;
            case 3
                ind = ind3;
        end
        
        % Simulate for 10V and store the result as the reference
        V = Data_PNG(ind(end));         % V -tThe voltage used
        
        [~,x_0] = ode45(@(t,x) Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk,Fzgrk,cN), tt, IC, opts);
        r10S = x_0(end,1)*1e3;      % m-->mm
        z10S = x_0(end,2)*1e3;      % m-->mm
        
        % Simulate for other voltages
        for kk = 1:length(ind)-1    %[V r z]
            V = Data_PNG(ind(kk),1);
            
            [~,x_0] = ode45(@(t,x) Robotics_Practicals_AL_ode(t,x,w_u,w_th,Amp_u,Amp_th,phi_u,phi_th,a,V,RRgrk,ZZgrk,Frgrk,Fzgrk,cN), tt, IC, opts);
            rS = x_0(end,1)*1e3;      % m-->mm
            zS = x_0(end,2)*1e3;      % m-->mm
                
            % Compute the differenc and relative error
            Sim_diff = sqrt( (r10S-rS)^2 + (z10S-zS)^2);
            Exp_diff = sqrt( (Data_PNG(ind(end),2)-Data_PNG(ind(kk),2))^2 + (Data_PNG(ind(end),3)-Data_PNG(ind(kk),3))^2);
            RErr(Sim_ind) = abs(Exp_diff - Sim_diff)/abs(Sim_diff);
                       
            fprintf('Particle %i, iter %i of %i, Err = %f\n',P_ind,kk,length(ind)-1,RErr(Sim_ind));
            Sim_ind = Sim_ind+1;
        end
    end
    RErr(isnan(RErr)) = [];     % Ignore unstable points

    Val = mean(RErr);       	% Compute the average
    figure(2); clf
        plot(RErr,'o')
        yline(Val,'--r')
        xlabel index
        ylabel Err
        title(sprintf('prescaler = %f, RErr = %f',prescaler,Val))

    drawnow; shg
end

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

