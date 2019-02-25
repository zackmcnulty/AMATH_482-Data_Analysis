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


%% Camera 1 case 1
plots = [0 0 0 1 0 0]; % which plots to show; called in get_xy_coords

close all; clc; 

video = vidFrames1_1;
xrange = [300,400];
yrange = [200,450];
var_scale = 1;
max_pixel_val = 240;

[x1_1, y1_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);

pause(100)

% Camera 2 case 1
close all; clc; 

video = vidFrames2_1;
xrange = [250, 350];
yrange = [100, 375];
var_scale = 1;
max_pixel_val = 240;

[x2_1, y2_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



% Camera 3 case 1
close all; clc; 

video = vidFrames3_1;
xrange = [250, 500];
yrange = [225, 325];
var_scale = 1;
max_pixel_val = 240;

[x3_1, y3_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);

%% Principal Component Analysis part 1
close all; clc;

rank_approx = 1;

% Some videos are behind the others in time; offset accounts for this
% by aligning the frames
% offset gives which from to start from for video 1,2, and 3 respectively.
offset = [11, 20, 11]; 
offset = offset - (min(offset) - 1);
pcs = 2; % number of principal components to plot
yrange = [100, 600];
A = my_pca(rank_approx, pcs, offset, yrange, x1_1, y1_1, x2_1, y2_1, x3_1, y3_1);


%% Part 2: Nosiy Case

clear all; close all; clc


load('camera_files/cam1_2.mat')
load('camera_files/cam2_2.mat')
load('camera_files/cam3_2.mat')


%% Camera 1 part 2
close all; clc; 

video = vidFrames1_2;
xrange = [300, 400];
yrange = [225, 400];
var_scale = 0.5;
max_pixel_val = 230;
plots = [0 0 0 0 0 0];

[x1_2, y1_2] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);


% Camera 2 part 2
close all; clc; 

video = vidFrames2_2;
xrange = [175, 450];
yrange = [50, 450];
var_scale = 0.5;
max_pixel_val = 240;
plots = [0 0 0 0 0 0];

[x2_2, y2_2] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



% Camera 3 part 2
close all; clc; 

video = vidFrames3_2;
xrange = [250, 500];
yrange = [225, 300];
var_scale = 0.8;
max_pixel_val = 230;
plots = [0 0 0 0 0 0];

[x3_2, y3_2] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




%% Principal Component Analysis part 2
close all; clc;

rank_approx = 3;
offset = [15, 1, 17];
offset = offset - (min(offset) - 1);
pcs = 4;
yrange = [100, 600];
A = my_pca(rank_approx, pcs, offset, yrange, x1_2, y1_2, x2_2, y2_2, x3_2, y3_2);



%% Part 3: Horizontal Displacement

clear all; close all; clc;

load('camera_files/cam1_3.mat')
load('camera_files/cam2_3.mat')
load('camera_files/cam3_3.mat')


%% Camera 1 part 3
close all; clc; 

video = vidFrames1_3;
xrange = [250, 400];
yrange = [200, 400];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 0 0];

[x1_3, y1_3] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




% Camera 2 part 3
close all; clc; 

video = vidFrames2_3;
xrange = [200, 400];
yrange = [175, 400];
var_scale = 1;
max_pixel_val = 240;
plots = [0 0 0 0 0 0];

[x2_3, y2_3] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




% Camera 3 part 3
close all; clc; 

video = vidFrames3_3;
xrange = [250, 450];
yrange = [175, 325];
var_scale = 1;
max_pixel_val = 245;
plots = [0 0 0 0 0 0];

[x3_3, y3_3] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);


%% Principal Component Analysis part 3

close all; clc;

rank_approx = 2;

% frame bucket moves to swinger's right
offset = [18 44 9];
offset = offset - (min(offset) - 1);
pcs = 3;
yrange = [100,600];
my_pca(rank_approx, pcs, offset, yrange, x1_3, y1_3, x2_3, y2_3, x3_3, y3_3);


%% Part 4: Horizontal Displacement AND Rotation
clear all; close all; clc;

load('camera_files/cam1_4.mat')
load('camera_files/cam2_4.mat')
load('camera_files/cam3_4.mat')

%% Camera 1 part 4

close all; clc; 

video = vidFrames1_4;
xrange = [300, 450];
yrange = [225, 400];
var_scale = 1;
max_pixel_val = 245;
plots = [0 0 0 0 0 0];

[x1_4, y1_4] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



% Camera 2 part 4
close all; clc; 

video = vidFrames2_4;
xrange = [210, 400];
yrange = [100, 350];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 0 0];

[x2_4, y2_4] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);


% Camera 3 part 4
close all; clc; 

video = vidFrames3_4;
xrange = [300, 500];
yrange = [175, 250];
var_scale = 0.7;
max_pixel_val = 235;
plots = [0 0 0 0 0 0];

[x3_4, y3_4] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Principal Component Analysis part 4
close all; clc;

rank_approx = 2;
offset = [11, 17, 9];
offset = offset - (min(offset) - 1);
pcs = 3;
yrange = [100,600];
my_pca(rank_approx, pcs, offset, yrange, x1_4, y1_4, x2_4, y2_4, x3_4, y3_4);
