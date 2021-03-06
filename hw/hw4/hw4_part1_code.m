% AMATH 482 HW 4

% Zachary McNulty

% NOTES: benefits of alignment is that we can drastically reduce the number of modes
% we have. SVD SUCKS with translated data so the uncropped images raise the supposed
% rank of the system becasue SVD cannot recognized the translated images as
% the same. Takes time for preprocessing.

% Convert of music into frequency space to avoid the issues with translation that interferes
% With SVD


%% Part 1: Yale Faces CROPPED

clear all; close all; clc;

% path to hw folder
% system('cd ~/Desktop/AMATH_482/hw/hw4')

% cropped vs uncropped shows how poor SVD does with translated images and
% the importance of pre-processing.
folders = dir('~/Desktop/AMATH_482/hw/hw4/input_files/CroppedYale/yale*');

% MY FACES!
F_ave = zeros(32256, 38); 
F = zeros(32256, 2415);

index = 1;
index2 = 1;
for j = 1:length(folders)
    path = strcat(folders(j).folder , '/', folders(j).name, '/*.pgm'); 
    files = dir(path);
    
    for k = 1:length(files)
        face = imread(strcat(files(k).folder , '/', files(k).name));
        
        % average all the images for each person
        % Since all the images have the same positional location of the
        % face this seems to help balance out the lighting in the
        % photograph.
        F(:, index2) = reshape(double(face), size(face, 1)*size(face,2), 1);
        index2 = index2 + 1;
        
        if k == 1
           ave_face = double(face); 
        else
            ave_face = ave_face + double(face);
        end    
    end
    ave_face = ave_face ./ length(files);
    
    % each row is a measurement of a given pixel across many trials (i.e.
    % photos of people)
    F_ave(:, index) = reshape(ave_face, size(ave_face, 1)*size(ave_face, 2), 1);
    index = index + 1;
end

% Demean the data
mean_ave = mean(F_ave, 2);
F_ave = F_ave - mean_ave;
mean_F = mean(F,2);
F = F - mean_F;

imheight = size(face,1);
imwidth = size(face,2);

%% Calculate the Singular Value Decomposition

%[u, s, v] = svd(F_ave, 'econ');
[u2, s2, v2] = svd(F, 'econ');


%% Interpretation of U, Sigma, V

% set which group of images we are working with, averaged or not
U = u;
S = s;
V = v;
X = F_ave;
m = mean_ave;

%%
% The columns of U represent an orthogonal basis for the codomain of our data matrix X.
% Not just any basis, however, the BEST basis for the space where each
% successive column captures as much of the variance in the data as
% possible.
% In this case, as the codomain of our matrix is the space in which column of
% the matrix exists, as each column is a face in our data set the codomain is 
% "face space" Thus, it follows that
% U is the "optimal" basis for face-space where the SVD's idea of what the important features
% that make up a face are will be those that capture the most variance.
% Furthermore, the columns of V then store the
% coordinates of each of our data measurements (i.e. each face) within this
% face space. What linear combination of these principal components
% (outputted by the code below).
num_faces = 5;
for col = 1:num_faces
   figure(2)
   
   subplot(1, num_faces, col)
   next_face = reshape(U(:, col), [imheight, imwidth]);
   pcolor(flipud(next_face)), shading interp, colormap(gray);
   title(strcat("Principal Component ", num2str(col)));
   set(gca,'xticklabel',[])
   set(gca,'yticklabel',[])
   set(gca, 'fontsize', 15)
end



%% Plot singular Values

figure(4)
plot(diag(S) ./ max(diag(S)), 'r.', 'markersize', 20)
%xticks(1:round(length(S)/10):length(S))
xticks(0:20:240)
xlim([0, 250])
ylim([0,1])
title('Normalized Singular Values')
ylabel('\sigma_j')
xlabel('index (j)')
set(gca, 'fontsize', 20)

%% Choosing rank r

% We will choose the r that captures at least 90% of the energy of the
% system

energy = 0;
total = sum(diag(S));
% how much energy we want our modes to capture.
% 75% does alright; 90% does very good.
threshold = 0.30; 
r = 0;
while energy < threshold
    r = r + 1;
    energy = energy + S(r,r)/total;
end


%% rank r approximation of face space

% low-rank approximation of our faces stored in F.
X_r = U(:, 1:rank) * S(1:rank, 1:rank) * (V(:, 1:rank)');

faces_to_plot = 10;

for j = 1:faces_to_plot
figure (5)
subplot(211)
imshow(uint8(reshape(X_r(:, j) + m, imheight, imwidth)));
title(strcat("Rank ", num2str(rank), " approximation"))
set(gca, 'fontsize', 15);
subplot(212)
imshow(uint8(reshape(X(:, j) + m, imheight, imwidth)));
title('Original Image')
set(gca, 'fontsize', 15);
end


%% Part 1: Yale Faces UNCROPPED

clear all; close all; clc;

% NOTE: in this case averaging might not work as well as the previous case
% because the photos are no longer aligned.

% cropped vs uncropped shows how poor SVD does with translated images and
% the importance of pre-processing.
folders = dir('~/Desktop/AMATH_482/hw/hw4/input_files/yalefaces_uncropped/subject*');

% MY FACES!
F_ave = zeros(77760, 15); 
F = zeros(77760, 165);

index = 1;
index2 = 1;
for j = 1:length(folders)
    path = strcat(folders(j).folder , '/', folders(j).name, '/subject*'); 
    files = dir(path);
    
    for k = 1:length(files)
        face = imread(strcat(files(k).folder , '/', files(k).name));
        
        % average all the images for each person
        % Since all the images have the same positional location of the
        % face this seems to help balance out the lighting in the
        % photograph.
        F(:, index2) = reshape(double(face), size(face, 1)*size(face,2), 1);
        index2 = index2 + 1;
        
        if k == 1
           ave_face = double(face); 
        else
           ave_face = ave_face + double(face);
        end    
    end
    ave_face = ave_face ./ length(files);
    
    % each row is a measurement of a given pixel across many trials (i.e.
    % photos of people)
    F_ave(:, index) = reshape(ave_face, size(ave_face, 1)*size(ave_face, 2), 1);
    index = index + 1;
end

% Demean data?
mean_ave = mean(F_ave, 2);
F_ave = F_ave - mean_ave;
mean_f = mean(F, 2);
F = F - mean(F);

imheight = size(face,1);
imwidth = size(face,2);

%% Calculate the Singular Value Decomposition

[u, s, v] = svd(F_ave, 'econ');
[u2, s2, v2] = svd(F, 'econ');


%% Interpretation of U, Sigma, V

% set which group of images we are working with, averaged or not
U = u2;
S = s2;
V = v2;
X = F;
m = mean_f;


%%
% The columns of U represent an orthogonal basis for the codomain of our data matrix X.
% Not just any basis, however, the BEST basis for the space where each
% successive column captures as much of the variance in the data as
% possible.
% In this case, as the codomain of our matrix is the space in which column of
% the matrix exists, as each column is a face in our data set the codomain is 
% "face space" Thus, it follows that
% U is the "optimal" basis for face-space where the SVD's idea of what the important features
% that make up a face are will be those that capture the most variance.
% Furthermore, the columns of V then store the
% coordinates of each of our data measurements (i.e. each face) within this
% face space. What linear combination of these principal components
% (outputted by the code below).

num_faces = 5;
for col = 1:num_faces
   figure(2)
   
   subplot(1, num_faces, col)
   next_face = reshape(U(:, col), [imheight, imwidth]);
   pcolor(flipud(next_face)), shading interp, colormap(gray);
   title(strcat("Principal Component ", num2str(col)));
   set(gca,'xticklabel',[])
   set(gca,'yticklabel',[])
   set(gca, 'fontsize', 15)
end




%% Plot singular Values

figure(4)
plot(diag(S) ./ max(diag(S)), 'r.', 'markersize', 20);
ylim([0,1])
title('Normalized Singular Values')
ylabel('\sigma_j')
xlabel('index (j)')
set(gca, 'fontsize', 20);
%% Choosing rank r

% We will choose the r that captures at least 90% of the energy of the
% system

energy = 0;
total = sum(diag(S));
% how much energy we want our modes to capture.
% 75% does alright; 90% does very good.
threshold = 0.9; 
r = 0;
while energy < threshold
    r = r + 1;
    energy = energy + S(r,r)/total;
end


%% rank r approximation of face space

rank = r;

% low-rank approximation of our faces stored in F.
X_r = U(:, 1:rank) * S(1:rank, 1:rank) * (V(:, 1:rank)') + m;

faces_to_plot = 10;

for j = 1:faces_to_plot
figure (5)
subplot(211)
imshow(uint8(reshape(X_r(:, j), imheight, imwidth)));
title('Low rank approximation')
subplot(212)
imshow(uint8(reshape(X(:, j) + m, imheight, imwidth)));
title("Original Image")
pause(1);
end







