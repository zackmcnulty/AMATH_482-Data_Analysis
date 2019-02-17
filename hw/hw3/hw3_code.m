%% HW 3: Principal Component Analysis (PCA)

% NOTE: convert uint8 to double using double() before processing!
% NOTE: each frame of video should only produce a single timepoint. Taken
%       the mean of the x and y values you find!
% get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots)

%% Part 1: Ideal Case


clear all; close all; clc;

% load data
load('camera_files/cam1_1.mat')
load('camera_files/cam2_1.mat')
load('camera_files/cam3_1.mat')


% Video loading taken from page 120 of class notes with slight
% modifications

%% Camera 1 case 1
close all; clc; 

video = vidFrames1_1;
xrange = [300,400];
yrange = [200,450];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 1 1];

[x1_1, y1_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Camera 2 case 1
close all; clc; 

video = vidFrames2_1;
xrange = [250, 350];
yrange = [100, 375];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 1 1];

[x2_1, y2_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Camera 3 case 1
close all; clc; 

video = vidFrames3_1;
xrange = [250, 500];
yrange = [225, 325];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 1 1];

[x3_1, y3_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);








%% Part 2: Nosiy Case

clear all; close all; clc



load('camera_files/cam1_2.mat')
load('camera_files/cam2_2.mat')
load('camera_files/cam3_2.mat')


%% Camera 1 part 2

%% Camera 2 part 2

%% Camera 3 part 2




%% Part 3: Horizontal Displacement

clear all; close all; clc;

load('camera_files/cam1_3.mat')
load('camera_files/cam2_3.mat')
load('camera_files/cam3_3.mat')


%% Camera 1 part 3

%% Camera 2 part 3

%% Camera 3 part 3




%% Part 4: Horizontal Displacement AND Rotation
clear all; close all; clc;

load('camera_files/cam1_4.mat')
load('camera_files/cam2_4.mat')
load('camera_files/cam3_4.mat')

%% Camera 1 part 4

%% Camera 2 part 4

%% Camera 3 part 4


