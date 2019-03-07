% AMATH 482 HW 4

% Zachary McNulty

% NOTES: benefits of alignment is that we can drastically reduce the number of modes
% we have. SVD SUCKS with translated data so the uncropped images raise the supposed
% rank of the system becasue SVD cannot recognized the translated images as
% the same. Takes time for preprocessing.


%% Part 1: Yale Faces CROPPED

clear all; close all; clc;

% path to hw folder
% system('cd ~/Desktop/AMATH_482/hw/hw4')

% cropped vs uncropped shows how poor SVD does with translated images and
% the importance of pre-processing.
folders = dir('~/Desktop/AMATH_482/hw/hw4/input_files/CroppedYale/yale*');

% MY FACES!
F = zeros(32256, 38); 
F2 = zeros(32256, 2415);

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
        F2(:, index2) = reshape(double(face), size(face, 1)*size(face,2), 1);
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
    F(:, index) = reshape(ave_face, size(ave_face, 1)*size(ave_face, 2), 1);
    index = index + 1;
end

% Demean the data
% F = F - mean(F,1);
% F2 = F2 - mean(F2, 1);

imheight = size(face,1);
imwidth = size(face,2);

%% Calculate the Singular Value Decomposition

[u, s, v] = svd(F, 'econ');
%[u2, s2, v2] = svd(F2, 'econ');


%% Interpretation of U, Sigma, V

% set which group of images we are working with, averaged or not
U = u;
S = s;
V = v;
X = F;

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

for col = 1:size(U, 2)
   figure(2)
   
   next_face = reshape(U(:, col), [imheight, imwidth]);
   pcolor(flipud(next_face)), shading interp, colormap(gray);

   pause(1);
end


%% Plot singular Values

figure(4)
plot(diag(S) ./ max(diag(S)), 'r.', 'markersize', 20)
xticks(1:round(length(S)/10):length(S))
ylim([0,1])
title('Normalized Singular Values')
ylabel('\sigma_j')
xlabel('index (j)')

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
X_r = U(:, 1:rank) * S(1:rank, 1:rank) * (V(:, 1:rank)');

faces_to_plot = 10;

for j = 1:faces_to_plot
figure (5)
subplot(211)
imshow(uint8(reshape(X_r(:, j), imheight, imwidth)));
subplot(212)
imshow(uint8(reshape(X(:, j), imheight, imwidth)));
pause(5);
end


%% Part 1: Yale Faces UNCROPPED

clear all; close all; clc;

% NOTE: in this case averaging might not work as well as the previous case
% because the photos are no longer aligned.

% cropped vs uncropped shows how poor SVD does with translated images and
% the importance of pre-processing.
folders = dir('~/Desktop/AMATH_482/hw/hw4/input_files/yalefaces_uncropped/subject*');

% MY FACES!
F = zeros(77760, 15); 
F2 = zeros(77760, 165);

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
        F2(:, index2) = reshape(double(face), size(face, 1)*size(face,2), 1);
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
    F(:, index) = reshape(ave_face, size(ave_face, 1)*size(ave_face, 2), 1);
    index = index + 1;
end

% Demean data?
% F = F - mean(F);
% F2 = F2 - mean(F2);

imheight = size(face,1);
imwidth = size(face,2);

%% Calculate the Singular Value Decomposition

[u, s, v] = svd(F, 'econ');
[u2, s2, v2] = svd(F2, 'econ');


%% Interpretation of U, Sigma, V

% set which group of images we are working with, averaged or not
U = u2;
S = s2;
V = v2;
X = F2;


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

for col = 1:size(U, 2)
   figure(2)
   
   next_face = reshape(U(:, col), [imheight, imwidth]);
   pcolor(flipud(next_face)), shading interp, colormap(gray);

   %pause(2);
end



%% Plot singular Values

figure(4)
plot(diag(S) ./ max(diag(S)), 'r.', 'markersize', 20)
xticks(1:round(length(S)/10):length(S))
ylim([0,1])
title('Normalized Singular Values')
ylabel('\sigma_j')
xlabel('index (j)')

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
X_r = U(:, 1:rank) * S(1:rank, 1:rank) * (V(:, 1:rank)');

faces_to_plot = 100;

for j = 1:faces_to_plot
figure (5)
subplot(211)
imshow(uint8(reshape(X_r(:, j), imheight, imwidth)));
title('Low rank approximation')
subplot(212)
imshow(uint8(reshape(X(:, j), imheight, imwidth)));
title("Original Image")
pause(1);
end

%%






%% Part 2: Music Classification


%% Test case 1: Different Bands from Different Genres


% Collect Data
clear all; close all; clc;

% path to hw folder
folders = dir('~/Desktop/AMATH_482/hw/hw4/input_files/music_files/part1/*');

% number of 5 second long samples to be taken from each song
% NOTE: takes actually twice this number of samples as half the samples are
% placed in the training set and the other half in the validation set.
samples_per_song = 5;
num_songs = 12; % total number of songs across all groups
sample_length = 5; % in seconds

% Fs is always 44100 across all the songs we downloaded based on
% the download procedure
Fs = 44100;


% Each row is a sample
% each data sample will be the frequency content of a given song
training_data = zeros(samples_per_song * num_songs, Fs*sample_length);
training_labels = zeros(samples_per_song * num_songs, 1);
validation_data = zeros(samples_per_song * num_songs, Fs*sample_length);
validation_labels = zeros(samples_per_song * num_songs, 1);

index = 1;
group_num = 1;
tic
for j = 1:length(folders)
    
    % skip all hidden folders within the directory
    if startsWith(folders(j).name, '.')
       continue 
    end
    
    path = strcat(folders(j).folder , '/', folders(j).name, '/*.mp3'); 
    files = dir(path);
    
    for k = 1:length(files)
        song_path = strcat(files(k).folder, "/", files(k).name);
        

        start = 30 * Fs; % skip first 30 seconds of song.
        finish = inf; % sample until the end of audio file
        
        % Fs is the sampling rate and Y is the amplitude at each point in
        % the recording
        [Y, Fs] = audioread(song_path, [start, finish]);
        
        % Y is a stereo measurement (includes measurements for both left
        % and right speakers) so we average these two to get a single
        % measurement
        
        Y = mean(Y, 2).';
        %p8 = audioplayer(Y,Fs); playblocking(p8);
       
        
        % randomly make training data set and cross validation data set
        for s = 1:samples_per_song
            
            % randomly choose a sample starting point
            sample_start = randi(length(Y) - sample_length*Fs, [2,1]);
            
            % from the given starting point, take a sample of
            % 'sample_length' seconds and store the corresponding label
            training_data(index, :) = fft(Y(sample_start(1):sample_start(1) + sample_length * Fs-1));
            training_labels(index) = group_num; % i.e. this song came from band j
            
            % take another sample to be used for validation
            validation_data(index, :) = fft(Y(sample_start(2):sample_start(2) + sample_length * Fs-1));
            validation_labels(index) = group_num;
            index = index + 1;
            
%             p8 = audioplayer(Y(sample_start(1):sample_start(1) + sample_length * Fs-1),Fs); playblocking(p8);
%             pause(3);
        end
        
    end
    
    group_num = group_num + 1;
    
end
toc




%% Find Dominant Modes for Each Group using SVD (i.e. each band)

% Extract the individual groups (i.e. bands in this case) from the training data
group1 = training_data(training_labels == 1);
group2 = training_data(training_labels == 2);
group3 = training_data(training_labels == 3);

% SVD the entire dataset to find principal components
% of "music space"
mean_td = mean(training_data, 2);
[u, s, v] = svd(training_data - mean_td, 'econ');

%% Plot singular values: which are relevant?

plot(diag(s) / max(diag(s)), 'r.', 'markersize', 20)
title('Normalized Singular Values')
ylabel('\sigma_j')
xlabel('index j')

% Find the modes that best separate the data by projecting each of our
% groups onto the individual components. Note that each column of V gives
% the coordinates of each data measurement (row in X = training_data) onto
% the corresponding principal component (column in U). So entry (i,j) of V
% gives the weighting/"importance" of principal component 




% Multiclass SVM
% https://www.mathworks.com/help/stats/classificationecoc-class.html
% let yoyoyo know if I get above 75% accuracy lol.