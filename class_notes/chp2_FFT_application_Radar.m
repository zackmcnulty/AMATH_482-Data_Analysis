

%% 2.2 FFT application in Radars

clear all; close all;
L=30; % Total time slot to transform
n=512; % number of Fourier modes 2^7; i.e. how many frequencies should I check???

% notice how we strip off the last point? Why?
% because we have periodic boundary conditions, so the last point will
% match the first one anyways.
t2=linspace(-L,L,n+1); t=t2(1:n); % time discretization 
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; % frequency components of FFT
u=sech(t); % ideal signal in the time domain 
figure(1), subplot(3,1,1), plot(t,u,'k'), hold on
xlim([-10,10])


noise=1;
ut=fft(u); % calculate the true frequency spectrum from ideal signal
utn=ut+noise*(randn(1,n)+i*randn(1,n)); % add gaussian distributed white noise to freq spectrum (i.e. random frequencies)
un=ifft(utn); % undo the fft to get the noisy signal.
figure(1), subplot(3,1,2), plot(t,abs(un),'k'), hold on
xlim([-10,10])



% increased noise;
noise=3;
ut=fft(u); % calculate the true frequency spectrum from ideal signal
utn=ut+noise*(randn(1,n)+i*randn(1,n)); % add gaussian distributed white noise to freq spectrum (i.e. random frequencies)
un=ifft(utn); % undo the fft to get the noisy signal.
figure(1), subplot(3,1,3), plot(t,abs(un),'k'), hold on
xlim([-10,10])


% Spectral filtering: filtering for freqs in a specific range
% When the noise gets too great, it is near impossible to determine true
% frequency? Filter for a specific range; 
% k0 = desired freq
% tau = selectivity (how close to k0 should freqs be? > tau --> more
% selective filter)

% generate noisy data in frequency domain
noise=10;
ut=fft(u); 
unt=ut+noise*(randn(1,n)+i*randn(1,n)); % add frequency noise
un=ifft(unt); % convert back to time series (w/ added freq noise)
subplot(2,1,1), plot(t,abs(un),'k') % plot time series of signal
axis([-30 30 0 2])

xlabel('time (t)'), ylabel('|u|')
subplot(2,1,2) 
plot(fftshift(k),abs(fftshift(unt))/max(abs(fftshift(unt))),'k') % plot freq domain of signal
axis([-25 25 0 1])
xlabel('wavenumber (k)')
ylabel('|ut|/max(|ut|)')

% apply a filter to the data, looking for freqs around wavenumber k0 = 0
% filter out freqs too far above or below this
filter=exp(-0.2*(k).^2); % tau = 0.2, k0 = 0
unft=filter.*unt;
unf=ifft(unft);

% figure 10 pg 39 of notes
figure
subplot(311)
plot(fftshift(k),abs(fftshift(unt))/max(abs(fftshift(unt))),'k') 
axis([-25 25 0 1])
title("unfiltered signal (frequency domain)")
xlabel("wavenumber")
ylabel("|u| / max(|u|)")


subplot(312)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'k') 
axis([-25 25 0 1])
title("frequency filtered signal (frequency domain)")
xlabel("wavenumber")
ylabel("|uf| / max(|uf|)")

subplot(313)
plot(t, abs(unf), 'k'), hold on;
plot(t, sech(t), 'r')
axis([-25 25 0 1])
legend({"filtered time signal", "ideal time signal (no noise)"})
title("time series data from filtered freqs")
xlabel('time(t)')
ylabel('|u|')



%% 2.3 FFT Application: Radar Detection and Averaging

clear all; close all; clc

% figure 12 pg 43 in class notes

figure
L=30; % total time slot
n=512; % Fourier modes
t2=linspace(-L,L,n+1); t=t2(1:n); 
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
ks=fftshift(k); 
noise=10;

labels=['(a)';'(b)';'(c)';'(d)']; 
realize=[1 2 5 100];

for jj=1:length(realize)
    u=sech(t); sum=zeros(1,n); ut=fft(u);
    
    % collect data for realize(jj) time frames and average the 
    % frequences collected across all of them
    for j=1:realize(jj)
        
        % add noise to freq domain
        utn(j,:)=ut+noise*(randn(1,n)+i*randn(1,n)); 
        sum=sum+utn(j,:); 
        dat(j,:)=abs(fftshift(utn(j,:)))/max(abs(utn(j,:))); 
        un(j,:)=ifft(utn(j,:));
    end
    
    % convert sum of frequencies to time course data and average it
    ave=abs(fftshift(sum))/realize(jj);
    subplot(4,1,jj) 
    plot(ks,ave/max(ave),'k')

    set(gca,'Fontsize',[15])
    axis([-20 20 0 1]) 
    text(-18,0.7,labels(jj,:),'Fontsize',[15]) 
    ylabel('|fft(u)|','Fontsize',[15])
end
hold on 

plot(ks,abs(fftshift(ut))/max(abs(ut)),'r:','Linewidth',[2]) 
set(gca,'Fontsize',[15])
xlabel('frequency (k)')


% figure 13, pg 44
figure(2)
waterfall(ks,1:8,dat(1:8,:)), colormap([0 0 0]), view(-15,80) 
set(gca,'Fontsize',[15],'Xlim',[-28 28],'Ylim',[1 8]) 
xlabel('frequency (k)'), ylabel('realization'),zlabel('|fft(u)|')




%% 2.4: Time-Frequency Analysis: Windowed Fourier Transforms

% no code!


%% 2.5 Time-Frequency Anaylsis and Wavelets!


