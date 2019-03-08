%% Part 2: Music Classification


%% Test case 1: Different Bands from Different Genres


% Collect Data
clear all; close all; clc;

% path to hw folder
folders = dir('~/Desktop/AMATH_482/hw/hw4/input_files/music_files/part1/*');

% number of 5 second long samples to be taken from each song
% NOTE: takes actually twice this number of samples as half the samples are
% placed in the training set and the other half in the validation set.
samples_per_song = 10;
num_songs = 21; % total number of songs across all groups
sample_length = 5; % in seconds
skip_period = 15; % skip the first "skip period" seconds of song to avoid aampling silence at beginning of song.
%sample_rate = 2; % take every other "sample_rate"th point.

% Fs is always 44100 across all the songs we downloaded based on
% the download procedure
Fs = 44100;


% Each row is a sample
% each data sample will be the frequency content of a given song
training_data = zeros(Fs*sample_length + 1, samples_per_song * num_songs);
%training_data = zeros(262152, samples_per_song * num_songs); % spectrogram
training_labels = zeros(samples_per_song * num_songs, 1);
validation_data = zeros(Fs*sample_length + 1, samples_per_song * num_songs);
%validation_data = zeros(262152, samples_per_song * num_songs); %spectro
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
        

        start = skip_period * Fs; % skip first "skip_period" seconds of song.
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
            % Take the fft of the sample to convert to frequency space
            
            training_data(:, index) = fft(Y(sample_start(1):sample_start(1) + sample_length * Fs)).';
            %spec = spectrogram(Y(sample_start(1):sample_start(1) + sample_length * Fs));
            %training_data(:, index) = reshape(spec, size(spec,1) * size(spec,2), 1);
            training_labels(index) = group_num; % i.e. this song came from band j
            
            % take another sample to be used for validation
            %validation_data(:, index) = fft(Y(sample_start(2):sample_start(2) + sample_length * Fs)).';
            
            %spec = spectrogram(Y(sample_start(2):sample_start(2) + sample_length * Fs));
            %validation_data(:, index) = reshape(spec, size(spec,1) * size(spec,2), 1);
            validation_labels(index) = group_num;
            index = index + 1;
        end
        
    end
    
    group_num = group_num + 1;
    
end
toc


%% Find Modes which are best at separating groups.

% Extract the individual groups (i.e. bands in this case) from the training data
group1 = training_data(:, training_labels == 1);
group2 = training_data(:, training_labels == 2);
group3 = training_data(:, training_labels == 3);

% SVD the entire dataset to find principal components
% of "music space"
mean_td = mean(training_data, 1);
[u, s, v] = svd(training_data - mean_td, 'econ'); % mean subtract

% Plot singular values: which are relevant?

singular_values = diag(s) / max(diag(s));
plot(singular_values, 'r.', 'markersize', 20)
title('Normalized Singular Values')
ylabel('\sigma_j')
xlabel('index j')

% Find the modes that best separate the data by projecting each of our
% groups onto the individual components. Note that each column of V gives
% the coordinates of each data measurement (row in X = training_data) onto
% the corresponding principal component (column in U). So entry (i,j) of V
% gives the weighting/"importance" of principal component 


%% Low rank approximation

% Calculate the low-rank approximation that captures at least 90% of the
% energy of the system.
energy_threshold = 0.99;
rank = find(cumsum(singular_values / sum(singular_values)) >= energy_threshold, 1);

X_r = u(:, 1:rank) * s(1:rank, 1:rank) * (v(:, 1:rank)') + mean_td;

%% Test Quality of reconstruction

% test = real(ifft(X_r(:, 1))); % all complex parts are zero and interfere with playing music
% p8 = audioplayer(test, Fs); 
% playblocking(p8);

%% Train a SVM model
clc;

% SVM Cannot handle complex numbers so we simply take the absolute value
% To convert the training data to amplitude (all real).

% In this case, as our data measurements are stored as columns, the columns of
% U store the fundamental structure/modes of our system (they form a basis
% for our domain : frequency / music space) while the columns of SV' give the
% coordinates of each data measurement (i.e. each song) within this
% frequency space. Thus, these coordinates define the underlying structure
% of our system and we can use them to classify.

% Specifically, we will use V* as our "coordinates" and project future
% data measurements onto the columns of U*S

vstar = v';
xtrain = [real(vstar(1:rank, :)) ; imag(vstar(1:rank, :))].';

%SVM model
Mdl_svm = fitcecoc(xtrain, training_labels);

% K nearest neighbors model
Mdl_knn = fitcknn(xtrain, training_labels);

% Naive Bayes Model
Mdl_nb = fitcnb(xtrain, training_labels);

% Decision Tree Classification model
Mdl_tree = fitctree(xtrain, training_labels);


% Computes the classification Error rate in classification for our training
% data set using the classification conditions defined by our model.
error_rate_training_svm = resubLoss(Mdl_svm)
error_rate_training_knn = resubLoss(Mdl_knn)
error_rate_training_nb = resubLoss(Mdl_nb)
error_rate_training_tree = resubLoss(Mdl_tree)

%% Test how good the model is at predicting labels


% since our model was trained on the coordinates given by V which are the coordinates
% of our data with respect to the basis US (for frequency space), we have to
% to get the coordinates of our validation data with respect to this
% basis: US * coordinates_in_US = data. Since U is unitary, this yields
% simply: coordinates_in_US = S^-1 U' * data
Sinv = diag(1 / diag(s));
coordinates_in_US = (u')*validation_data;
xtest = [real(coordinates_in_US(1:rank, :)) ; imag(coordinates_in_US(1:rank, :))].';

% we only used the first r = rank coordinates to determine our
% classification so we again do that.
predictions_svm = predict(Mdl_svm, xtest);
predictions_knn = predict(Mdl_knn, xtest);
predictions_nb = predict(Mdl_nb, xtest);
predictions_tree = predict(Mdl_tree, xtest);
accuracy_svm = sum(predictions_svm == validation_labels) / length(validation_labels)
accuracy_knn = sum(predictions_knn == validation_labels) / length(validation_labels)
accuracy_nb = sum(predictions_nb == validation_labels) / length(validation_labels)
accuracy_tree = sum(predictions_tree == validation_labels) / length(validation_labels)






