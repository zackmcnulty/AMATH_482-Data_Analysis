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
plots = [0 0 0 0 0 0];

[x1_1, y1_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Camera 2 case 1
close all; clc; 

video = vidFrames2_1;
xrange = [250, 350];
yrange = [100, 375];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 0 0];

[x2_1, y2_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Camera 3 case 1
close all; clc; 

video = vidFrames3_1;
xrange = [250, 500];
yrange = [225, 325];
var_scale = 1;
max_pixel_val = 240;
plots = [0 0 0 0 0 0];

[x3_1, y3_1] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);

%% Principal Component Analysis part 1
close all; clc;

rank_approx = 1;
A = my_pca(rank_approx, x1_1, y1_1, x2_1, y2_1, x3_1, y3_1);

figure(2)
subplot(221)
plot(y1_1)
subplot(222)
plot(x1_1)
subplot(223)
plot(A(2,:))
subplot(224)
plot(A(1,:))



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
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 0 0 0 0];

[x1_2, y1_2] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);


%% Camera 2 part 2
close all; clc; 

video = vidFrames2_2;
xrange = [175, 450];
yrange = [50, 450];
var_scale = 0;
max_pixel_val = 250;
plots = [0 0 0 0 0 0];

[x2_2, y2_2] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




%% Camera 3 part 2
close all; clc; 

video = vidFrames3_2;
xrange = [250, 500];
yrange = [225, 300];
var_scale = 0.1;
max_pixel_val = 245;
plots = [0 0 0 0 0 0];

[x3_2, y3_2] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




%% Principal Component Analysis part 2
close all; clc;

my_pca(x1_2, y1_2, x2_2, y2_2, x3_2, y3_2);





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
plots = [0 0 1 0 1 1];

[x1_3, y1_3] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




%% Camera 2 part 3
close all; clc; 

video = vidFrames2_3;
xrange = [200, 400];
yrange = [175, 400];
var_scale = 1;
max_pixel_val = 240;
plots = [0 0 1 0 1 1];

[x2_3, y2_3] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);




%% Camera 3 part 3
close all; clc; 

video = vidFrames3_3;
xrange = [250, 450];
yrange = [175, 325];
var_scale = 1;
max_pixel_val = 245;
plots = [0 0 1 0 1 1];

[x3_3, y3_3] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);


%% Principal Component Analysis part 3

close all; clc;

my_pca(x1_3, y1_3, x2_3, y2_3, x3_3, y3_3);


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
plots = [0 0 1 1 1 1];

[x1_4, y1_4] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Camera 2 part 4
close all; clc; 

video = vidFrames2_4;
xrange = [210, 400];
yrange = [100, 350];
var_scale = 1;
max_pixel_val = 250;
plots = [0 0 1 1 1 1];

[x2_4, y2_4] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Camera 3 part 4
close all; clc; 

video = vidFrames3_4;
xrange = [300, 500];
yrange = [175, 250];
var_scale = 0.7;
max_pixel_val = 235;
plots = [0 0 1 1 1 1];

[x3_4, y3_4] = get_xy_coords(video, xrange, yrange, var_scale, max_pixel_val, plots);



%% Principal Component Analysis part 2
close all; clc;

my_pca(x1_4, y1_4, x2_4, y2_4, x3_4, y3_4);
