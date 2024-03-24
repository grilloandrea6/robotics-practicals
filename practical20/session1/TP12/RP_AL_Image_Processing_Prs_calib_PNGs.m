%% Robotics Practicals - Image processing stored png files
% Use image processing to identify the instantanous position of 
% the particle and the position and oriantation of the levitator
%
% Use 2D cross-correlation + subpixel resolution to identify the position 
% of the particle in space and markers on the levitator
%
% You are asked to load the calibration file of the camera and the
% estimated radii of the markers from the rotation axis. Then you are asked
% to load a single *.avi video to be processed.
%
% This script is used to estimate:
% 1) The position and oriantation of the levitator
% 2) The position of the particle
% 2) The particles' radii

clear; clc; close all
SAVE = 1;                       % 1 - save data, 0 - don't save data
Filenew = 'Calib_PNG_data.mat'; % The script saves this file
magnification = 2e2;    % magnification coeficient
% Load camera parameters
[File, path] = uigetfile('*.mat','Select the camera calibration results','MultiSelect','off');
load([path File]);

% Load dvice radii
[File, path] = uigetfile('*.mat','Load the device radii','MultiSelect','off');
load([path File]);

% Load images
[Files, path] = uigetfile('*.png','Select all the calibration images','MultiSelect','ON');

%% Set the matrices
Nf = length(Files(:));                          % No. of images
Pos_Mat = nan*zeros(Nf,2);                      % Prepare a matrix to store the position of the particle
Pos_Mat_TG = Pos_Mat;                           % Prepare a matrix to store the position of the top marker
Pos_Mat_BG = Pos_Mat;                           % Prepare a matrix to store the position of the bottom marker
Particle_Rad_mmV = zeros(Nf,1);                 % Prepare a vector to store the position of the bottom marker
O_angle_rad = zeros(Nf,1);                      % store the device instantanous angle (rad)
V_V = zeros(Nf,1);                              % store the applied voltage (V)
P_index = zeros(Nf,1);                          % store the particle index

%% Process the first image
A0 = imread([path Files{1}]);                   % load an image and undistort it
A0 = undistortImage(A0,cameraParams); 
A_lim = [0 max(A0(:))];

figure(1);clf  
    imshow(double(A0),A_lim, 'InitialMagnification', magnification); hold all

    title('Choose the Top region (hollow marker) to estimate the oriantation (paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixTO = TEMP(2):TEMP(2)+TEMP(4); iyTO = TEMP(1):TEMP(1)+TEMP(3);

    title('Choose the Bottom region (full marker) to estimate the oriantation (paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixBO = TEMP(2):TEMP(2)+TEMP(4); iyBO = TEMP(1):TEMP(1)+TEMP(3);

    title('Choose the general region to track the particle (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixGP = TEMP(2):TEMP(2)+TEMP(4); iyGP = TEMP(1):TEMP(1)+TEMP(3);

figure(1);clf  
    imshow(double(A0),A_lim, 'InitialMagnification', magnification);

    title('Choose a small region around the hollow marker (Top) (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixTOc = TEMP(2):TEMP(2)+TEMP(4); iyTOc = TEMP(1):TEMP(1)+TEMP(3);

    title('Choose a small region around the full marker (Bottom) (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixBOc = TEMP(2):TEMP(2)+TEMP(4); iyBOc = TEMP(1):TEMP(1)+TEMP(3);

    title('Choose a small region around the particle (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixSP = TEMP(2):TEMP(2)+TEMP(4); iySP = TEMP(1):TEMP(1)+TEMP(3);

OriantationT = double(A0(ixTO,iyTO));           % Oriantation Top region
OriantationB = double(A0(ixBO,iyBO));           % Oriantation Bottom region

OriantationTc = double(A0(ixTOc,iyTOc));        % Oriantation Top circle region
OriantationBc = double(A0(ixBOc,iyBOc));        % Oriantation Bottom circle region
SmallP = double(A0(ixSP,iySP));                 % Small particle tracking region

%% Process the images 
for kk = 1:Nf
    
    % Request applied voltage and particle index
    prompt = {'Voltage value', 'Particle number' };
    dlgtitle = Files{kk};
    dims = [1 100]; definput = {'10','1'}; answer = inputdlg(prompt,dlgtitle,dims,definput);
    V_V(kk) = str2double(answer{1});        % Real Frame rate
    P_index(kk) = str2double(answer{2});    % Hz - Excitation frequency
    
    % Coose the region to estimate the location of the top hollow marker
    % Coose the region to estimate the location of the bottom full marker
    % Choose the general region to estimate the particle location
    % Choose small regions aroud the markers and particle
    
    A0 = imread([path Files{kk}]);                      % load an image and undistort it
    A0 = undistortImage(A0,cameraParams); 

    if kk>=2 &&  P_index(kk)~=P_index(kk)               % If  different particle is considered update the region being tracked
        figure(1);clf 
        title('Choose a small region around the particle (Paused)')
        ROI = drawrectangle; pause;
        TEMP = round(ROI.Position); ixSP = TEMP(2):TEMP(2)+TEMP(4); iySP = TEMP(1):TEMP(1)+TEMP(3);
        SmallP = double(A0(ixSP,iySP));                 % Small particle tracking region
    end
    
    % Estimate the particle's radius
    GeneralP = double(A0(ixGP,iyGP));                   % General particle tracking region
    % Estimate the particle radius 
    figure(2); clf 
    imshow(GeneralP,A_lim, 'InitialMagnification', magnification*100);
    drawnow
    title('Draw a circle around the particle to estimate its radius (paused)','FontSize',9);
    ROI = drawcircle; pause;
    Particle_Rad_mmV(kk) = ROI.Radius*Pix2mm;
    % Particle position                                  
    COR2D = normxcorr2(SmallP,GeneralP);                                    % 2D correlation
    Pos_pix = peakfit2d_ad(COR2D);                                          % find max corolation with sub-pixel interpolation
   
    Pos_Mat(kk,:) = [Pos_pix(2)-0.5*size(SmallP,1),0.5*size(SmallP,2)-Pos_pix(1)] + [iyGP(1) -ixGP(1)];     % Save the position in the original frame   
      
    % Device position and oriantation 
    OriantationT = double(A0(ixTO,iyTO));
    OriantationB = double(A0(ixBO,iyBO));
    COR2DT = normxcorr2(OriantationTc,OriantationT);                        % 2D correlation
    COR2DB = normxcorr2(OriantationBc,OriantationB);                        % 2D correlation
   
    Pos_pixT = peakfit2d_ad(COR2DT);                                        % find max with sub-pixel interpolation 
    Pos_pixB = peakfit2d_ad(COR2DB);                                        % find max with sub-pixel interpolation 

    Pos_pixTG = Pos_pixT + [ixTO(1)-0.5*size(OriantationTc,2) iyTO(1)-0.5*size(OriantationTc,1)];   % in the original frame 
    Pos_pixBG = Pos_pixB + [ixBO(1)-0.5*size(OriantationBc,2) iyBO(1)-0.5*size(OriantationBc,1)];   % in the original frame
    O_angle_rad(kk) = atan2( (-Pos_pixTG(2)+Pos_pixBG(2)),(-Pos_pixTG(1)+Pos_pixBG(1)) );
    Pos_Mat_TG(kk,:) = [Pos_pixTG(2),-Pos_pixTG(1)];
    Pos_Mat_BG(kk,:) = [Pos_pixBG(2),-Pos_pixBG(1)];
    
    radiusT_px = radiusT/Pix2mm;
    radiusB_px = radiusB/Pix2mm;
    Norm = sqrt(sum((Pos_Mat_BG(kk,:)-Pos_Mat_TG(kk,:)).^2,2));  
    CT = Pos_Mat_TG(kk,:)+radiusT_px*(Pos_Mat_BG(kk,:)-Pos_Mat_TG(kk,:))./Norm;     % The center as measured from the top marker
    CB = Pos_Mat_BG(kk,:)+radiusB_px*(Pos_Mat_TG(kk,:)-Pos_Mat_BG(kk,:))./Norm;     % The center as measured from the bottom marker
    C = 0.5*(CT+CB);                                                                % Average of the center location
    
    figure(2);
    title([int2str(kk) ' [' int2str(Nf) ']' ]), shg, drawnow

    drawnow; shg;

end
%% Post process the data
Pos_Mat_TG_mm = Pix2mm*Pos_Mat_TG;      % Px-->mm
Pos_Mat_BG_mm = Pix2mm*Pos_Mat_BG;      % Px-->mm
Pos_Mat_mm = Pos_Mat * Pix2mm;          % Px-->mm
% 1) Get the origin position (u,v) and track the particle in (r,z) with respect to the origin
Norm = sqrt(sum((Pos_Mat_BG_mm-Pos_Mat_TG_mm).^2,2));               % Norm of the vector connecting the top and bottom markers

CT_v = Pos_Mat_TG_mm+radiusT*(Pos_Mat_BG_mm-Pos_Mat_TG_mm)./Norm;   % The center as measured from the top marker
CB_v = Pos_Mat_BG_mm+radiusB*(Pos_Mat_TG_mm-Pos_Mat_BG_mm)./Norm;   % The center as measured from the bottom marker
C_v = 0.5*(CT_v+CB_v);                                              % Average of the center location

r = (Pos_Mat_mm(:,1)-C_v(:,1)).*cos(O_angle_rad) + ...
    (Pos_Mat_mm(:,2)-C_v(:,2)).*sin(O_angle_rad);

z = (C_v(:,1)-Pos_Mat_mm(:,1)).*sin(O_angle_rad) + ...
    (Pos_Mat_mm(:,2)-C_v(:,2)).*cos(O_angle_rad);

%% Sort according to particle index
[P_index,sort_i] = sort(P_index);
V_V = V_V(sort_i);
C_v(:,:) = C_v(sort_i,:);
r(:) = r(sort_i,:);
z(:) = z(sort_i,:);
O_angle_rad = O_angle_rad(sort_i);
Particle_Rad_mmV = Particle_Rad_mmV(sort_i);

%% Sort acording to voltage
for kk = 1:length(unique(P_index))
    [V_V(P_index==kk),sort_i] = sort(V_V(P_index==kk));
    Temp = C_v(P_index==kk,:);
    C_v(P_index==kk,:) = Temp(sort_i,:);
    Temp = r(P_index==kk,:);
    r(P_index==kk) = Temp(sort_i);
    Temp = z(P_index==kk,:);
    z(P_index==kk) = Temp(sort_i);
    Temp = O_angle_rad(P_index==kk);
    O_angle_rad(P_index==kk) = Temp(sort_i);
    Temp = Particle_Rad_mmV(P_index==kk);
    Particle_Rad_mmV(P_index==kk) = Temp(sort_i);
end

%% plot
ind1 = find(P_index==1);
ind2 = find(P_index==2);
ind3 = find(P_index==3);
aN(3,1) = mean(Particle_Rad_mmV(ind3));
aN(2,1) = mean(Particle_Rad_mmV(ind2));
aN(1,1) = mean(Particle_Rad_mmV(ind1));
O_angle_radN(3,1) = mean(O_angle_rad(ind3));
O_angle_radN(2,1) = mean(O_angle_rad(ind2));
O_angle_radN(1,1) = mean(O_angle_rad(ind1));

figure(3); clf
    yyaxis left
    h(1,1) = plot(V_V(ind1),z(ind1),'--o'); hold all
    h(2,1) = plot(V_V(ind2),z(ind2),'--x'); 
    xlabel('V (V)')
    ylabel('z (mm)')
    yyaxis right
    h(1,1) = plot(V_V(ind1),r(ind1),'--o'); 
    h(2,1) = plot(V_V(ind2),r(ind2),'--x'); 
    ylabel('r (mm)')
    

%% Save
if SAVE==1  
    save([path  Filenew],'aN','O_angle_radN','r','z','V_V','ind1','ind2','ind3');
end


%% Functions
function POS = peakfit2d_ad(Z)
%% interpolate the location of sub-pizel maximum
% whithin a matrix Z representing the cross correlation
% COR2D = normxcorr2(TRACK_region,GEN_region);    % 2D correlation
% Pos_pix = peakfit2d_ib(COR2D)

% peak approximation using 2D polynomial fit within 9 point neighboring points hip
% find global maximum and extract 9 neighbor points
[~,IND] = max(Z(:));
[ii,jj] = ind2sub(size(Z),IND); 
K = Z(ii-1:ii+1,jj-1:jj+1); % 9 points around maximal pixel
%%      z = c+b'*[x;y]+1/2*[x y]*B*[x;y]  - model to fit
x = -1:1; 
y = -1:1;
[X,Y] = meshgrid(x,y);          % create grid
A = [0.5*X(:).^2 Y(:).*X(:) X(:) 0*X(:)+1 0.5*Y(:).^2 Y(:)]; % model
aa = A\K(:);                    %Solve (Least squares)
B = [aa(1) aa(2); aa(2) aa(5) ];% B matrix in model
b = -[aa(3);aa(6)];             % b vecto
h = B\b;                        % solution for dz/dx=0, dz/dy=0

POS =[h(2)+ii h(1)+jj];  % interpolated location (y,x)
end

function J = Linear_sine_fit(f,tt,signal)
w = 2*pi*f;
% Model
A = [cos(w*tt) sin(w*tt) cos(2*w*tt) sin(2*w*tt) cos(3*w*tt) sin(3*w*tt)];% 1+0*tt tt];
a = A\signal;   % Solve LS problem

J = rms(A*a-signal); % compute the residual
end

function [center,radius,residuals] = circlefit(x,y)
% Fit a circle to data using the least squares approach.
%
% Fits the equation of a circle in Cartesian coordinates to a set of xy
% data points by solving the overdetermined system of normal equations, i.e.
% x^2 + y^2 + a*x + b*y + d = 0
%
% Input arguments:
% x,y
%    Cartesian coordinates of noisy data points
%
% Output arguments:
% center:
%    coordinates of the least-squares fit circle center
% radius:
%    least-squares fit cirlce radius
% residuals:
%    residuals in the radial direction
%
% Examples:
% [center,radius,residuals] = circlefit(X)
% [center,radius,residuals] = circlefit(x,y,z);

% Copyright 2010 Levente Hunyadi

narginchk(1,2);
switch nargin  % n x 2 matrix
    case 1
        n = size(x,1);
        validateattributes(x, {'numeric'}, {'2d','real','size',[n,2]});
        y = x(:,2);
        x = x(:,1);
    otherwise  % two x,y vectors
        validateattributes(x, {'numeric'}, {'real','vector'});
        validateattributes(y, {'numeric'}, {'real','vector'});
        n = numel(x);
        x = x(:);  % force into columns
        y = y(:);
        validateattributes(x, {'numeric'}, {'size',[n,1]});
        validateattributes(y, {'numeric'}, {'size',[n,1]});
end

% solve linear system of normal equations
A = [x, y, ones(size(x))];
b = -(x.^2 + y.^2);
a = A \ b;

% return center coordinates and sphere radius
center = -a(1:2)./2;
radius = realsqrt(sum(center.^2)-a(3));

if nargout > 2
	% calculate residuals
   residuals = radius - sqrt(sum(bsxfun(@minus,[x y],center.').^2,2));
end
end