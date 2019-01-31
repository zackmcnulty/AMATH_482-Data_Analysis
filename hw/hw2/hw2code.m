% AMATH 482 HW 2: Gabor Transforms

close all; clear all; clc

load handel
v = y'/2;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

%p8 = audioplayer(v,Fs);
%playblocking(p8);

%%
close all; clc;

v = y'/2;
v = v(1:end-1); % remove last point so we have even number of points
L = 9; % tspan(end); length of our signal in time (or space) domain


tspan = (1:length(v))/Fs;
n = length(v); % ALWAYS has to be even

% rescale [-pi pi] domain fft assumes to fit our domain [-L, L]
k = 2*pi / (2*L) .* [0:(n/2 - 1) -n/2:-1];
ks = fftshift(k);  

tstep = 0.5; % how much to slide the window after each measurement to new measurement
power = 2; % power of the gaussian 
a = 1; %width parameter of our gabor filter; actually width is prop to 1/a

% list of values to center filter at; collect frequencies for each
% one of the windows centered at this point for spectrogram
tau = 0:tstep:tspan(end); 

spectogram_data = zeros(n, length(tau));

figure(1)

for j = 1:length(tau)
    % Gaussian Window
    %gabor_filter = exp(-a*(tspan - tau(j)).^power);
    
    % Mexican Hat Window
    %gabor_filter = (2 ./ (sqrt(3*a)*pi^0.25)) .* (1 - ((tspan - tau(j)) ./ a).^2).*exp(-1.*(tspan - tau(j)).^2 ./ (2*a^2));
    
    % Shannon (Step-function) Window
    gabor_filter = abs(tspan - tau(j)) <= 1/(2*a);
    
    % plot our original signal
    subplot(411)
    plot((1:length(v))/Fs,v);
    ylim([-0.5,0.5])
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Signal of Interest, v(n)');
    
    % Plot our gabor filter
    subplot(412)
    plot(tspan, gabor_filter);
    ylim([-2, 3]) % mexican hat wavelet
    %ylim([0,1.1]) %gaussian limits
   
    xlim([0,9])
    xlabel('Time [sec]');
    title('Gabor Filter')
    ylabel('Filter Magnitude');
    
    % plot our filtered signal
    v_filtered = gabor_filter .* v;
    
    subplot(413)
    plot(tspan, v_filtered);
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Filtered Signal in Spatial Domain')
    ylim([-0.5,0.5])
    
    
    % plot our filtered signals frequency domain
    vft = fft(v_filtered);
    
    subplot(414)
    plot(ks, abs(fftshift(vft)) ./ max(abs(vft)));
    xlabel('Time [sec]');
    ylabel('|uft| / max(|uft|)');
    title('Filtered Signal in Frequency Domain')
    ylim([0, 1.2])
    
    
    drawnow
 
    spectogram_data(:, j) = (abs(fftshift(vft)) ./ max(abs(vft)) ).';
 
end

   % make spectogram plot
   figure(2)
   pcolor(tau, ks, spectogram_data), shading interp;
   ylabel('Frequency (wavenumber)')
   xlabel('Time [sec]')


   
   
   
   
   
   
   
   
   
   
   
%% Part 2: Piano
clear all; close all; clc;

 figure(1)
 tr_piano=16;  % record time in seconds
 m1 = audioread('music1.wav').';
 Fs1=length(m1)/tr_piano;
 plot((1:length(m1))/Fs1,m1);
 xlabel('Time [sec]'); ylabel('Amplitude');
 title('Mary had a little lamb (piano)');  drawnow
 %p8 = audioplayer(m1,Fs1); playblocking(p8);
 
 
 tspan1 = (1:length(m1))/Fs1;
 L1 = tr_piano;
 n1 = length(m1);
 k1 = 2*pi / (2*L1) .* [0:(n1/2 - 1) -n1/2:-1];
 ks1 = fftshift(k1);
 tstep1 = 0.5; % how much to slide the window after each measurement to new measurement
 power1 = 2; % power of the gaussian 
 a1 = 10; %width parameter of our gabor filter; higher a = lower width

% list of values to center filter/gabor window at
tau1 = 0:tstep1:tspan1(end); 

% stores the max wavenumber for each gabor window (centered at tau(j))
max_k1 = zeros(1, length(tau1));
 
for j = 1:length(tau1)
    gabor_filter1 = exp(-a1*(tspan1 - tau1(j)).^power1); % Gaussian
    
    % mexican hat
    %gabor_filter = (2 ./ (sqrt(3*a1)*pi^0.25)) .* (1 - ((tspan1 - tau1(j)) ./ a1).^2).*exp(-1.*(tspan1 - tau1(j)).^2 ./ (2*a1^2));
    
    % plot our original signal
    subplot(411)
    plot((1:length(m1))/Fs1,m1);
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Mary Had a little lamb (piano)');
    
    % Plot our gabor filter
    subplot(412)
    plot(tspan1, gabor_filter1);
    %ylim([-2, 3]) % mexican hat wavelet
    ylim([0,1.1]) %gaussian limits
   
    xlim([0,L1])
    xlabel('Time [sec]');
    title('Gabor Filter')
    ylabel('Filter Magnitude');
    
    % plot our filtered signal
    m1_filtered = gabor_filter1 .* m1;
    
    subplot(413)
    plot(tspan1, m1_filtered);
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Filtered Signal in Spatial Domain')
    ylim([-0.7,0.7]);
    xlim([0, L1])
    
    
    % plot our filtered signals frequency domain
    m1ft = fft(m1_filtered);
 
    % filter overtones?
    
    subplot(414)
    plot(ks1, abs(fftshift(m1ft)) ./ max(abs(m1ft)));
    xlabel('Wavenumber (k)');
    ylabel('|uft| / max(|uft|)');
    title('Filtered Signal in Frequency Domain')
    xlim([-20000,20000])
    
    % outside "max" call is to filter out the negative wavenumber
    max_k1(j) = max(k1(abs(m1ft) == max(abs(m1ft))));
   
    drawnow
  
end



% convert wavenumber to hertz! 
max_freqs1 = max_k1 ./ (2*pi);

%%

notes = ["A", "A#", "B", "C", "C#","D", "D#", "E", "F", "F#", "G", "G#"];
fund_freqs = [27.5, 29.135, 30.863, 32.703, 34.648,36.708, 38.891, 41.203, 43.654,46.249, 48.99, 51.913];

all_notes = [];
all_freqs = [];

for j = 1:8
    for note = notes
        all_notes = [all_notes  strcat(note, num2str(j))];
    end
    for freq = fund_freqs
        % freqs double each time you move up scale
        all_freqs = [all_freqs freq*(2^(j-1))];
    end
end

for freq = max_freqs1
    [min_val, index] = min(abs(freq - all_freqs));
    
    all_notes(index)
    freq
    all_freqs(index)
end


% find fundamental note whose frequency has the lowest remainder (i.e. this
% note is a multiple (overtone?) of that fundamental one

% for freq = max_freqs
%    [min_val, index] = min(abs(freq - round(freq./fund_freqs).*fund_freqs));
%    
%    notes(index)
%    fund_freqs(index)
% end


% SPECTROGRAM SHOWING JUST BLOCKS OF THE MOST COMMON FREQ!
% Expect EDC notes only?


%% Part 2: Recorder
clear all; close all; clc;


 figure(2)
 tr_rec=14;  % record time in seconds
 m2=audioread('music2.wav').'; 
 Fs2=length(m2)/tr_rec;
 plot((1:length(m2))/Fs2,m2);
 xlabel('Time [sec]'); ylabel('Amplitude');
 title('Mary had a little lamb (recorder)');
 %p8 = audioplayer(m2,Fs2); playblocking(p8);
 
 tspan2 = (1:length(m2))/Fs2;
 L2 = tr_rec;
 n2 = length(m2);
 k2 = 2*pi / (2*L2) .* [0:(n2/2 - 1) -n2/2:-1];
 ks2 = fftshift(k2);
 tstep2 = 0.5; % how much to slide the window after each measurement to new measurement
 power2 = 2; % power of the gaussian 
 a2 = 10; %width parameter of our gabor filter

% list of values to center filter at; collect frequencies for each
% one of the windows centered at this point for spectrogram
tau2 = 0:tstep2:tspan2(end); 

% stores the max wavenumber for each gabor window (centered at tau(j))
max_k2 = zeros(1, length(tau2));
 
for j = 1:length(tau2)
    gabor_filter2 = exp(-a2*(tspan2 - tau2(j)).^power2); % Gaussian
    
    % mexican hat
    %gabor_filter = (2 ./ (sqrt(3*a2)*pi^0.25)) .* (1 - ((tspan2 - tau2(j)) ./ a2).^2).*exp(-1.*(tspan2 - tau2(j)).^2 ./ (2*a2^2));
    
    % plot our original signal
    subplot(411)
    plot((1:length(m2))/Fs2,m2);
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Mary Had a little lamb (recorder)');
    
    % Plot our gabor filter
    subplot(412)
    plot(tspan2, gabor_filter2);
    %ylim([-2, 3]) % mexican hat wavelet
    ylim([0,1.1]) %gaussian limits
   
    xlim([0,L2])
    xlabel('Time [sec]');
    title('Gabor Filter')
    ylabel('Filter Magnitude');
    
    % plot our filtered signal
    m2_filtered = gabor_filter2 .* m2;
    
    subplot(413)
    plot(tspan2, m2_filtered);
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Filtered Signal in Spatial Domain')
    ylim([-0.7,0.7]);
    xlim([0, L2])
    
    
    % plot our filtered signals frequency domain
    m2ft = fft(m2_filtered);
 
    % filter overtones?
    
    subplot(414)
    plot(ks2, abs(fftshift(m2ft)) ./ max(abs(m2ft)));
    xlabel('Wavenumber (k)');
    ylabel('|uft| / max(|uft|)');
    title('Filtered Signal in Frequency Domain')
    xlim([-20000,20000])
    
    % outside "max" call is to filter out the negative wavenumber
    max_k2(j) = max(k2(abs(m2ft) == max(abs(m2ft))));
   
    drawnow
  
end



% convert wavenumber to hertz! 
max_freqs2 = max_k2 ./ (2*pi);


notes = ["A", "A#", "B", "C", "C#","D", "D#", "E", "F", "F#", "G", "G#"];
fund_freqs = [27.5, 29.135, 30.863, 32.703, 34.648,36.708, 38.891, 41.203, 43.654,46.249, 48.99, 51.913];

all_notes = [];
all_freqs = [];

for j = 1:8
    for note = notes
        all_notes = [all_notes  strcat(note, num2str(j))];
    end
    for freq = fund_freqs
        % freqs double each time you move up scale
        all_freqs = [all_freqs freq*(2^(j-1))];
    end
end
 
for freq = max_freqs2
    [min_val, index] = min(abs(freq - all_freqs));
    
    all_notes(index)
    freq
    all_freqs(index)
end
