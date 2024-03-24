%% Robotics Practicals - Save calibration parameters to file
% This code is used to save the calibration parameters obtained in the 
% camera calibration app as well as to estimate the px to mm ratio in the
% acoustic transducer plane.

% Initiate the script
clc; clear; close all

%Input measured checkerboard square size
squareSizeMm = 4.583;

%Initialize cameraParams variable and open app
cameraParams = [];
cameraCalibrator

%Wait for parameters to be created
while isempty(cameraParams); pause(0.1); end

%Image in ball plane
[Files, path] = uigetfile('*.bmp','Select coplanar image','MultiSelect','Off');
SampleImage = imread([path Files]);

%undistort with calibration params
[newIm, ~] = undistortImage(SampleImage,cameraParams);

%Find coordinates of checkerboard corners
[imagePointsPlanar,~] = detectCheckerboardPoints(newIm);

%Calculate distance between neighboring points in px
squareSizePx = sqrt((imagePointsPlanar(1,1)-imagePointsPlanar(2,1))^2+(imagePointsPlanar(1,2)-imagePointsPlanar(2,2))^2);

%Get ratio mm/px
Pix2mm = squareSizeMm/squareSizePx; %mm/px

save('calibrationData.mat','Pix2mm','cameraParams');