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
%sample_rate = 2; % take every other "sample_rate"th point.

% Fs is always 44100 across all the songs we downloaded based on
% the download procedure
Fs = 44100;


% Each row is a sample
% each data sample will be the frequency content of a given song
training_data = zeros(samples_per_song * num_songs, Fs*sample_length + 1);
training_labels = zeros(samples_per_song * num_songs, 1);
validation_data = zeros(samples_per_song * num_songs, Fs*sample_length + 1);
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
            % Take the fft of the sample to convert to frequency space
            training_data(index, :) = fft(Y(sample_start(1):sample_start(1) + sample_length * Fs));
            training_labels(index) = group_num; % i.e. this song came from band j
            
            % take another sample to be used for validation
            validation_data(index, :) = fft(Y(sample_start(2):sample_start(2) + sample_length * Fs));
            validation_labels(index) = group_num;
            index = index + 1;
            
            % p8 = audioplayer(Y(sample_start(1):sample_start(1) + sample_length * Fs-1),Fs); playblocking(p8);
            % pause(3);
        end
        
    end
    
    group_num = group_num + 1;
    
end
toc


%% Find Modes which are best at separating groups.

% Extract the individual groups (i.e. bands in this case) from the training data
group1 = training_data(training_labels == 1, :);
group2 = training_data(training_labels == 2, :);
group3 = training_data(training_labels == 3, :);

% rearrange training data to be ordered by group
% training_data(:, :, :) = [group1; group2; group3];

% SVD the entire dataset to find principal components
% of "music space"
mean_td = mean(training_data, 2);
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
energy_threshold = 0.90;
rank = find(cumsum(singular_values / sum(singular_values)) >= energy_threshold, 1);
X_r = u(:, 1:rank) * s(1:rank, 1:rank) * (v(:, 1:rank)');

%% Test Quality of reconstruction

test = real(ifft(X_r(2,:))); % all complex parts are zero and interfere with playing music
p8 = audioplayer(test, Fs); 
playblocking(p8);

%% Train a SVM model
clc;

% SVM Cannot handle complex numbers so we simply take the absolute value
% To convert the training data to amplitude (all real).

% In this case, as our data measurements are stored as rows, the columns of
% V store the fundamental structure/modes of our system (they form a basis
% for our domain : frequency / music space) while the columns of U*Sigma give the
% coordinates of each data measurement (i.e. each song) within this
% frequency space. Thus, these coordinates define the underlying structure
% of our system and we can use them to classify.

Mdl = fitcecoc(abs(u(1:rank, :).'), training_labels);
% Mdl = fitcecoc(abs(training_data), training_labels);
% Mdl = fitcecoc(abs(X_r), training_labels);


% Computes the classification Error rate in classification for our training
% data set using the classification conditions defined by our model.
error_rate_training = resubLoss(Mdl)

%% Test how good the model is at predicting labels

correct_predictions = 0;
for j = 1:length(validation_labels)
    % since our model was trained on the coordinates given by U which are the coordinates
    % of our data with respect to the basis V (for frequency space), we have to
    % to get the coordinates of our validation data with respect to this
    % basis: V * coordinates_in_V = data. Since V is unitary, this yields
    % simply: coordinates_in_V = V' * data
    
    next_song = validation_data(j, :).';
    coordinates_in_V = v'* next_song;
    
    % we only used the first r = rank coordinates to determine our
    % classification so we again do that.
    prediction = predict(Mdl,abs(coordinates_in_V(1:rank)).'); 
    actual = validation_labels(j);
    if prediction == actual
       correct_predictions = correct_predictions + 1 ;
    end
end

accuracy = correct_predictions / length(validation_labels)


% Multiclass SVM
% https://www.mathworks.com/help/stats/classificationecoc-class.html
% let yoyoyo know if I get above 75% accuracy lol.