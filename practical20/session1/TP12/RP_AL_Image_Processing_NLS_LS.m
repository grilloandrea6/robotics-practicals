%% Robotics Practicals - Controlling the position of an acoustically levitated particle
% Use image processing tools to identify the instantanous position of 
% the particle and the position and oriantation of the levitator

% Use cross correlation + subpixel resolution to identify the position of 
% the sphere in space and markers on the levitator

clear; clc; close all
SAVE = 0;               % 1 - save data, 0 - don't save data
magnification = 2e2;    % magnification coeficient
% Load camera parameters
%[File, path] = uigetfile('*.mat','Select the camera calibration results','MultiSelect','off');
%load([path File]);
load("calibrationData.mat")
load("Device_Radii.mat")
% Load dvice radii
%[File, path] = uigetfile('*.mat','Load the device radii','MultiSelect','off');
%load([path File]);

% Load a video file
%[File, path] = uigetfile('*.avi','Select a single video','MultiSelect','off');
% path = '/home/andrea'
% File = 

% Request prescalar value
% prompt = {'The real frame rate (FPS)', 'The exciatin frequency (Hz)' };
% dlgtitle = 'To remove the dialog box comment lines 15-20 & uncoment lines 22-24';
% dims = [1 100]; definput = {'50','1'}; answer = inputdlg(prompt,dlgtitle,dims,definput);
% Fs = str2double(answer{1});     % Real Frame rate
% Fex = str2double(answer{2});    % Hz - Excitation frequency

Fs = 50;                % Real Frame rate
Fex = 2;              % Hz - Excitation frequency

% load the video
s = '/home/andrea/2hz.avi';
v = VideoReader(s,'CurrentTime',0); % open video file

NpW = v.Width;                                  % No. of pixel w 
NpH = v.Height;                                 % No. of pixel H 
Nf = ceil(v.FrameRate*v.Duration)-1;            % No. of frames
Pos_Mat_TG = zeros(Nf,2);                       % Prepare a matrix to captur the position of the top marker
Pos_Mat_BG = zeros(Nf,2);                       % Prepare a matrix to captur the position of the bottom marker
A0 = undistortImage(readFrame(v),cameraParams); % Capture an image and undistort it
A0 = rgb2gray(A0);

A_lim = [0 max(A0(:))];

O_angle_rad = zeros(Nf,1);                      % store the device instantanous angle (rad)
%% Average the video to estimate the relevant regions for image processing
v = VideoReader(s,'CurrentTime',0);             % open video file
Atemp = 0;
while hasFrame(v)
    % Particle position
    AA = undistortImage(readFrame(v),cameraParams); % Capture an image and undistort it
    AA = rgb2gray(AA);
    Atemp = Atemp+(AA)/Nf;
end 
%% Process the first frame
% Coose the region to estimate the oriantation
% Choose the region to estimate the displacement
% Choose the general area to capture the particle displacement
% Choose small region aroud the particle
%figure(2);clf  
%    imshow(double(A0),A_lim, 'InitialMagnification', magnification); hold all
figure(1);clf  
    imshow(Atemp,double(A_lim)*0.2, 'InitialMagnification', magnification); hold all
 
    title('Choose the Top region (hollow marker) to estimate the oriantation (paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixTO = TEMP(2):TEMP(2)+TEMP(4); iyTO = TEMP(1):TEMP(1)+TEMP(3);
    
    title('Choose the Bottom region (full marker) to estimate the oriantation (paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixBO = TEMP(2):TEMP(2)+TEMP(4); iyBO = TEMP(1):TEMP(1)+TEMP(3);

figure(1);clf  
    imshow(double(A0),A_lim, 'InitialMagnification', magnification);
    
    title('Choose a small region around hollow marker (Top) (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixTOc = TEMP(2):TEMP(2)+TEMP(4); iyTOc = TEMP(1):TEMP(1)+TEMP(3);
    
    title('Choose a small region around full marker (Bottom) (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixBOc = TEMP(2):TEMP(2)+TEMP(4); iyBOc = TEMP(1):TEMP(1)+TEMP(3);

OriantationT = double(A0(ixTO,iyTO));           % Oriantation Top region
OriantationB = double(A0(ixBO,iyBO));           % Oriantation Bottom region

OriantationTc = double(A0(ixTOc,iyTOc));        % Oriantation Top circle region
OriantationBc = double(A0(ixBOc,iyBOc));        % Oriantation Bottom circle region

%% Process the Video
% Estimate the position of the particle, the device position and its angle for each frame
v = VideoReader(s,'CurrentTime',0);             % Open the video file
i_frame = 1;                                    % frame counter
while hasFrame(v)
    % Particle position
    Atemp = undistortImage(readFrame(v),cameraParams);                      % Read the next frame and undistort it
    Atemp = rgb2gray(Atemp);
    title([int2str(i_frame) ' [' int2str(Nf) ']' ]), shg, drawnow
   
    % Device position and oriantation 
    OriantationT = double(Atemp(ixTO,iyTO));
    OriantationB = double(Atemp(ixBO,iyBO));
    COR2DT = normxcorr2(OriantationTc,OriantationT);                        % 2D correlation
    COR2DB = normxcorr2(OriantationBc,OriantationB);                        % 2D correlation
    %%
    Pos_pixT = peakfit2d_ad(COR2DT);                                        % find max with sub-pixel interpolation 
    Pos_pixB = peakfit2d_ad(COR2DB);                                        % find max with sub-pixel interpolation 
    %%
    Pos_pixTG = Pos_pixT + [ixTO(1)-0.5*size(OriantationTc,2) iyTO(1)-0.5*size(OriantationTc,1)];   % in the original frame 
    Pos_pixBG = Pos_pixB + [ixBO(1)-0.5*size(OriantationBc,2) iyBO(1)-0.5*size(OriantationBc,1)];   % in the original frame
    O_angle_rad(i_frame) = atan2( (-Pos_pixTG(2)+Pos_pixBG(2)),(-Pos_pixTG(1)+Pos_pixBG(1)) );
    Pos_Mat_TG(i_frame,:) = [Pos_pixTG(2),-Pos_pixTG(1)];
    Pos_Mat_BG(i_frame,:) = [Pos_pixBG(2),-Pos_pixBG(1)];
    
%     figure(3);clf
%         imshow(OriantationT,A_lim, 'InitialMagnification', magnification); hold on
%         plot(Pos_pixT(2)-0.5*size(OriantationTc,1),Pos_pixT(1)-0.5*size(OriantationTc,2),'ro')
%         
%     figure(4);clf
%         imshow(OriantationB,A_lim, 'InitialMagnification', magnification); hold on
%         plot(Pos_pixB(2)-0.5*size(OriantationBc,1),Pos_pixB(1)-0.5*size(OriantationBc,2),'ro')

    
    i_frame = i_frame+1;  % increase frame counter
%     drawnow; shg;
%     pause(0.1)
end
timeV = 1/Fs * (0:i_frame-2)';
%% Post process the data
Pos_Mat_TG_mm = Pix2mm*Pos_Mat_TG;
Pos_Mat_BG_mm = Pix2mm*Pos_Mat_BG;
% 1) Fit harmonic signal to the origin (u,v)<->(x,y)
% 2) Fit harmonic signal to the instantanous angle (theta)
% 3) Track the particle in (r,z) with respect to the origin
%% 1) Fit harmonic signal to the origin (u,v)<->(x,y)
Norm = sqrt(sum((Pos_Mat_BG_mm-Pos_Mat_TG_mm).^2,2));               % Norm of the vector connecting the top and bottom markers

CT_v = Pos_Mat_TG_mm+radiusT*(Pos_Mat_BG_mm-Pos_Mat_TG_mm)./Norm;   % The center as measured from the top marker
CB_v = Pos_Mat_BG_mm+radiusB*(Pos_Mat_TG_mm-Pos_Mat_BG_mm)./Norm;   % The center as measured from the bottom marker
C_v = 0.5*(CT_v+CB_v);                                              % Average of the center location

t_ss = timeV(timeV>=0);
t_ss = t_ss - t_ss(1);
TSCALE = t_ss(end);
tau = t_ss/TSCALE;
U_pos_mm_ss = C_v(timeV>=0,1);   % For u - X coordinate
V_pos_mm_ss = C_v(timeV>=0,2);   % For V - Y coordinate

% U<->X coordinate
% Estimate the frequency using NL-LS
F_est = fminbnd(@(ff) Linear_sine_fit(ff,tau,detrend(U_pos_mm_ss)/max(abs(detrend(U_pos_mm_ss)))), 0.7*Fex*TSCALE, 1.3*Fex*TSCALE);
Fex_ident_u = F_est/TSCALE;
% Solve LS problem for input and output
w_est = 2*pi*Fex_ident_u;
Amat_uv = [cos(w_est*t_ss) sin(w_est*t_ss) cos(2*w_est*t_ss) sin(2*w_est*t_ss) cos(3*w_est*t_ss) sin(3*w_est*t_ss) 1+0*t_ss];
coefs_u = Amat_uv\U_pos_mm_ss;                  % Solve LS problem for the input
DC_u = coefs_u(end);
[Amp_u, IM] = max([sqrt(coefs_u(1)^2+coefs_u(2)^2) sqrt(coefs_u(3)^2+coefs_u(4)^2) sqrt(coefs_u(5)^2+coefs_u(6)^2)]);
Fex2_ident_u = Fex_ident_u*IM;

% V<->Y coordinate -  solve LS using the same frequency
coefs_v = Amat_uv\V_pos_mm_ss;                  % Solve LS problem for the input
DC_v = coefs_v(end);
Amp_v = max([sqrt(coefs_v(1)^2+coefs_v(2)^2) sqrt(coefs_v(3)^2+coefs_v(4)^2) sqrt(coefs_v(5)^2+coefs_v(6)^2)]);

UV_pos_mm_shift = C_v-[DC_u DC_v];


%% 2) Fit harmonic signal to the instantanous angle (theta)
O_angle_rad_ss = O_angle_rad(timeV>=0);
    
% Estimate the frequency using NL-LS
F_est = fminbnd(@(ff) Linear_sine_fit(ff,tau,detrend(O_angle_rad_ss)/max(abs(detrend(O_angle_rad_ss)))), 0.7*Fex*TSCALE, 1.3*Fex*TSCALE);
Fex_ident_th = F_est/TSCALE;
% Solve LS problem for input and output
w_est = 2*pi*Fex_ident_th;
Amat_th = [cos(w_est*t_ss) sin(w_est*t_ss) cos(2*w_est*t_ss) sin(2*w_est*t_ss) cos(3*w_est*t_ss) sin(3*w_est*t_ss) 1+0*t_ss];
    
coefs_th = Amat_th\O_angle_rad_ss;                  % Solve LS problem for the input
DC_th = coefs_th(end);
[Amp_th, IM] = max([sqrt(coefs_th(1)^2+coefs_th(2)^2) sqrt(coefs_th(3)^2+coefs_th(4)^2) sqrt(coefs_th(5)^2+coefs_th(6)^2)]);
Fex1_ident_th = Fex_ident_th*IM;
figure(4);clf
subplot(2,1,1)
plot(timeV,[UV_pos_mm_shift(:,1) Amat_uv*coefs_u-DC_u, U_pos_mm_ss-(Amat_uv*coefs_u)],'--o')
    ylabel('x (mm)')
    legend('Measured','Fit','Error')
    title(sprintf('Fitted frequency: %.3f Hz, Amplitude: %.3f mm, DC: %.3f mm',Fex2_ident_u,Amp_u,DC_u))
%subplot(2,1,2)
% plot(timeV,[UV_pos_mm_shift(:,2) Amat_uv*coefs_v-DC_v, V_pos_mm_ss-(Amat_uv*coefs_v)],'--o')
%     xlabel('t (s)')
%     ylabel('y (mm)')
%     title(sprintf('Amplitude: %.3f mm, DC: %.3f mm',Amp_v,DC_v))
%     legend('Measured','Fit','Error')

% figure(5);clf
subplot(2,1,2)
plot(timeV,[O_angle_rad Amat_th*coefs_th O_angle_rad-Amat_th*coefs_th]/pi*180,'--o'); hold all
yline(DC_th/pi*180,'r')
    xlabel('t (s)')
    ylabel('\theta (\circ)')
    legend('Measured','Fit','Error','DC offset')
    title(sprintf('Fitted frequency: %.3f Hz, DC offset: %.3f, Amplitude: %.3f ',Fex1_ident_th,DC_th/pi*180,Amp_th/pi*180))
figure(5); clf

if SAVE==1
    Filenew = strrep(File,'avi','mat');
    save([path  Filenew],'timeV','Fex_ident_u','U_pos_mm_ss','V_pos_mm_ss','Amp_v','Amp_u','DC_u','DC_v','Fex_ident_th','O_angle_rad','Amp_th','DC_th');
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