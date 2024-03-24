%% Robotics Practicals - Image processing recorded *avi videos
% Use image processing to identify the instantanous position of and 
% oriantation of the levitator
%
% Use 2D cross-correlation + subpixel resolution to identify the position 
% of the markers on the levitator in space 
%
% You are asked to load the calibration file of the camera. Then you are 
% asked to load a single *.mp4 video to be processed.
%
% This script is used to estimate:
% 1) The position and oriantation of the levitator
% 2) The radii of the marker' center relative to the rotation axis.

% Inititiate the script
clear; clc; close all
SAVE = 1;               % 1 - save data, 0 - don't save data
magnification = 2e2;    % magnification coeficient
% Load camera parameters
[File, path] = uigetfile('*.mat','Select the camera calibration results','MultiSelect','off');
load([path File]);

% Load a video file
[File, path] = uigetfile('*.mp4','Select a single video','MultiSelect','off');
s = [path File];
v = VideoReader(s,'CurrentTime',0); % open the video file

Nf = ceil(v.FrameRate*v.Duration)-1;            % No. of frames
Pos_Mat_TG = zeros(Nf,2);                       % Prepare a matrix to captur the position of the top marker
Pos_Mat_BG = zeros(Nf,2);                       % Prepare a matrix to captur the position of the bottom marker
O_angle_rad = zeros(Nf,1);                      % Prepare a vector to store the device instantanous angle (rad)

A0 = undistortImage(readFrame(v),cameraParams); % Capture an image and undistort it
A0 = rgb2gray(A0);
A_lim = [0 max(A0(:))];

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
% Coose the region to estimate the location of the top hollow marker
% Coose the region to estimate the location of the bottom full marker
% Choose small regions aroud the markers
figure(2);clf  
    imshow(double(A0),A_lim, 'InitialMagnification', magnification); hold all
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
    
    title('Choose a small region around the hollow marker (Top) (Paused)')
    ROI = drawrectangle; pause;
    TEMP = round(ROI.Position); ixTOc = TEMP(2):TEMP(2)+TEMP(4); iyTOc = TEMP(1):TEMP(1)+TEMP(3);
    
    title('Choose a small region around the full marker (Bottom) (Paused)')
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
    Atemp = undistortImage(readFrame(v),cameraParams); %  % Read the next frame and undistort it
    Atemp = rgb2gray(Atemp);
                                      
    title([int2str(i_frame) ' [' int2str(Nf) ']' ]), shg, drawnow
      
    % Device position and oriantation 
    OriantationT = double(Atemp(ixTO,iyTO));
    OriantationB = double(Atemp(ixBO,iyBO));
    COR2DT = normxcorr2(OriantationTc,OriantationT);                        % 2D correlation
    COR2DB = normxcorr2(OriantationBc,OriantationB);                        % 2D correlation

    Pos_pixT = peakfit2d_ad(COR2DT);                                        % find max with sub-pixel interpolation 
    Pos_pixB = peakfit2d_ad(COR2DB);                                        % find max with sub-pixel interpolation 

    Pos_pixTG = Pos_pixT + [ixTO(1)-0.5*size(OriantationTc,2) iyTO(1)-0.5*size(OriantationTc,1)];   % in the original frame 
    Pos_pixBG = Pos_pixB + [ixBO(1)-0.5*size(OriantationBc,2) iyBO(1)-0.5*size(OriantationBc,1)];   % in the original frame
    O_angle_rad(i_frame) = atan2( (-Pos_pixTG(2)+Pos_pixBG(2)),(-Pos_pixTG(1)+Pos_pixBG(1)) );
    Pos_Mat_TG(i_frame,:) = [Pos_pixTG(2),-Pos_pixTG(1)];
    Pos_Mat_BG(i_frame,:) = [Pos_pixBG(2),-Pos_pixBG(1)];
    
    i_frame = i_frame+1;                                                    % increase frame counter

end
%% Post process the data
Pos_Mat_TG_mm = Pix2mm*Pos_Mat_TG;      % Px-->mm
Pos_Mat_BG_mm = Pix2mm*Pos_Mat_BG;      % Px-->mm

% 1) Fit circles to the trajectories of the hollow and full markers
COL = lines(3);

[centerT,radiusT,~] = circlefit(Pos_Mat_TG_mm(:,1),Pos_Mat_TG_mm(:,2));
[centerB,radiusB,~] = circlefit(Pos_Mat_BG_mm(:,1),Pos_Mat_BG_mm(:,2));

Norm = sqrt(sum((Pos_Mat_BG_mm-Pos_Mat_TG_mm).^2,2));               % Norm of the vector connecting the top and bottom markers

CT_v = Pos_Mat_TG_mm+radiusT*(Pos_Mat_BG_mm-Pos_Mat_TG_mm)./Norm;   % The center as measured from the top marker
CB_v = Pos_Mat_BG_mm+radiusB*(Pos_Mat_TG_mm-Pos_Mat_BG_mm)./Norm;   % The center as measured from the bottom marker
C_v = 0.5*(CT_v+CB_v);                                              % Average of the center location

thvP = linspace(0,2*pi,1e3)';
xxT = radiusT.*cos(thvP)+centerT(1);
yyT = radiusT.*sin(thvP)+centerT(2);

xxB = radiusT.*cos(thvP)+centerB(1);
yyB = radiusT.*sin(thvP)+centerB(2);

figure(5);clf
h(1,1) = plot(Pos_Mat_TG_mm(:,1),Pos_Mat_TG_mm(:,2),'o','Color',COL(1,:)); hold all
h(2,1) = plot(xxT,yyT,'Color',COL(1,:)); 
h(3,1) = plot(centerT(1),centerT(2),'o','MarkerFaceColor',COL(1,:),'MarkerEdgeColor',COL(1,:)); 
h(4,1) = plot(Pos_Mat_BG_mm(:,1),Pos_Mat_BG_mm(:,2),'s','Color',COL(2,:)); 
h(5,1) = plot(xxB,yyB,'Color',COL(2,:));  
h(6,1) = plot(centerB(1),centerB(2),'s','MarkerFaceColor',COL(2,:),'MarkerEdgeColor',COL(2,:)); 
h(7,1) = plot(C_v(:,1),C_v(:,2),'x','MarkerFaceColor',COL(3,:),'MarkerEdgeColor',COL(3,:)); 

    axis equal
    grid minor
    hl = legend('Top marker - measurments','Fitted circle','Fitted center','Bottom marker - measurments','Fitted circle','Fitted center','Device center');
    set(hl,'Location','EastOutside')

%% 2) Plot the instantanous angle (theta)
figure(6);clf
plot(unwrap(O_angle_rad)/pi*180,'--o'); hold all
    xlabel('frame (#)')
    ylabel('\theta(\circ)')
    axis tight

%% 3) Plot the instantanous position
figure(7);clf
subplot(2,1,1)
plot(C_v(:,1),'--o'); hold all
    xlabel('frame (#)')
    ylabel('u (mm)')
    axis tight

subplot(2,1,2)
plot(C_v(:,2),'--o'); hold all
    xlabel('frame (#)')
    ylabel('v (mm)')
    axis tight

if SAVE==1
    Filenew = strrep(File,'mp4','mat');
    save([path  Filenew],'O_angle_rad','C_v');
    save([path  'Device_Radii.mat'],'radiusT','radiusB');
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