% AMATH 482 HW 5: Dynamic Mode Decomposition
% Zachary McNulty

% DMD for dynamic data?

%% Loading the video clips
clear all; close all; clc;
%video = VideoReader('input_files/car.MOV');
%video = VideoReader('input_files/newtcradle.mp4');
%video = VideoReader('input_files/mythbusters.mp4');
%video = VideoReader('input_files/watermelon.mp4');
video = VideoReader('input_files/bird2.mp4');

flip_contrast = 1; % flip contrast of pixels for analysis

%get(video) gives info on video file

frame_skip = 2; % take every "frame_rate"th frame
resolution_reduction = 0.2; % reduce resolution by this factor.
imheight = video.Height*resolution_reduction;
imwidth = video.Width*resolution_reduction;
num_frames = ceil(video.Duration * video.FrameRate / frame_skip);
dt = frame_skip / video.FrameRate; % used later
% t = 0: dt :video.Duration;
dt = 1;
t = 1:num_frames;

v = zeros(num_frames, imheight, imwidth); % stores frames as rows of 2D matrices
X = zeros(imheight*imwidth, num_frames); % stores frames as 1D columns


frame = 1;
index = 1;
while hasFrame(video)
   next_frame = readFrame(video);
   if mod(index, frame_skip) == 0
        x = imresize(next_frame, resolution_reduction);
        %imshow(x);
        if flip_contrast
            v(frame, :, :) = imcomplement(rgb2gray(x)); % imcomplement
        else
            v(frame, :, :) = rgb2gray(x); % dont imcomplement
        end
        X(:, frame) = reshape(v(frame,:,:), [imheight*imwidth, 1]);
        frame = frame + 1;
   end 
   index = index + 1;
end

% % demean data!
% mean_data = mean(mean(X));
% X = X - mean_data;

%% Make X and X' matrices
clc; close all;

%%%%%% body of DMD %%%%%%%%%%
X1 = X(:,1:end-1); X2 = X(:,2:end);

[U2,Sigma2,V2] = svd(X1, 'econ');

threshold = 0.90;
r = find(cumsum(diag(Sigma2) ./ sum(diag(Sigma2))) > threshold  ,1);
%r = size(Sigma2, 1); % full rank
figure(1)
subplot(121)
plot(diag(Sigma2)./max(diag(Sigma2)), 'r.', 'markersize', 20);
title('Normalized Singular Values')
xlabel('index ( j )')
ylabel('\sigma_j')
yticks(0:0.1:1)
set(gca, 'fontsize', 20);
subplot(122);
plot(log(diag(Sigma2)./max(diag(Sigma2))), 'r.', 'markersize', 20);
title('Normalized Singular Values')
xlabel('index ( j )')
ylabel('Log(\sigma_j )')
yticks(-8:1:0)
set(gca, 'fontsize', 20);


U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
Atilde = U'*X2*V/Sigma;
[W,D] = eig(Atilde);
Phi=X2*V/Sigma*W;

mu=diag(D);
omega=log(mu)/dt;

y0 = Phi\X(:, 1);  % pseudo-inverse initial conditions


[min_omega, min_index] = min(abs(omega));
foreground_modes_indices = find(abs(omega) > min_omega);
% omega(foreground_modes_indices)

%%
clc; close all;

X_lowrank = y0(min_index).*Phi(:, min_index).*exp(omega(min_index).*t);
X_sparse = X - abs(X_lowrank);

R = X_sparse .* (X_sparse < 0); % places all negative entries in R

% R = zeros(size(X_sparse)); % But why?

X_lowrank2 = abs(X_lowrank) + R; % R stores moving object? Cancel out from background and add to foreground
X_sparse2 = X_sparse - R; 

% Recreate the Foreground and Background

clc; close all;

figure(2)
set(gcf, 'Position',  [100, 100, 1000, 1000])
for frame = 1:size(X_lowrank2, 2)
    if flip_contrast
       lowrank = imcomplement(reshape(X_lowrank2(:, frame), [imheight, imwidth]));
       sparse = imcomplement(reshape(X_sparse2(:, frame), [imheight, imwidth] ));
       rj = imcomplement(reshape(R(:, frame), [imheight, imwidth]));
    else
        
       lowrank = reshape(X_lowrank2(:, frame), [imheight, imwidth]);
       sparse = reshape(X_sparse2(:, frame), [imheight, imwidth] );
       rj = reshape(R(:, frame), [imheight, imwidth]);
    end

   subplot(221)
   %imshow(uint8(lowrank));
   pcolor(flipud(lowrank)), shading interp, colormap(gray);
   title('Background (X\_lowrank)')
   set(gca, 'fontsize', 20);
   
   subplot(222)
   %imshow(uint8(sparse));
   pcolor(flipud(sparse)), shading interp, colormap(gray);
   title('Moving Subject (X\_sparse)')
   set(gca, 'fontsize', 20);
    
   subplot(223)
   %imshow(uint8(sparse + lowrank + rj));
   pcolor(flipud(sparse + lowrank)), shading interp, colormap(gray);
   title('Sum of both')
   set(gca, 'fontsize', 20);
   
   subplot(224)
   %imshow(uint8(abs(rj)));
   pcolor(flipud(rj)), shading interp, colormap(gray);
   title('R matrix')
   set(gca, 'fontsize', 20);
   
   drawnow;
end







