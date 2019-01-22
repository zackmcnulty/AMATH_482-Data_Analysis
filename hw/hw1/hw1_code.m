% HW 1: Saving my god damn dog

clear all; close all; clc;
load Testdata

L=15; % spatial domain 
n=64; % Fourier modes


x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

total_t = zeros(n,n,n);

for j=1:20 % we have 20 samples of data: time data
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   close all, isosurface(X,Y,Z,abs(Un))
   axis([-20 20 -20 20 -20 20]), grid on, drawnow
   title("Unfiltered Signal")
   pause(1)
end

%% Averaging of the spectrum

% we have 20 samples of data: time data
% use these different realizations of the object/data to
% average out zero mean noise.

% sum up all frequency data
for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   total_t = total_t + fftn(Un);
end

% find max value so we can normalize
flat = reshape(total_t, 1, n^3);
M = max(abs(flat));

% normalize the data and shift in the frequency domain so
% it is centered around zero.
ave = abs(fftshift(total_t)) ./ M;

% since we normalized ave, the max must be one.
% find frequency index where max occurs
ind2 = find(ave == 1);

% central frequencies in x,y, and z directions respectively.
xc = Kx(ind2);
yc = Ky(ind2);
zc = Kz(ind2);


%% Filtering of the spectrum

tau = 0.3; % filtering bandwidth

% create 3D gaussian filter. If you are far from central x,y, OR z
% frequency, your signal will be filtered out. This is in shifted
% frequency space however.
filter = exp(-1*tau * ((Kx - xc).^2 +  (Ky - yc).^2 + (Kz - zc).^2));

% unshift filter as it is currently centered at zero.
filter_s = ifftshift(filter);


% apply this filter to frequency domain at each timestep
% to find marbles location at each step
marble_coords = zeros(20, 3);

for j=1:20 % we have 20 samples of data: time data
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   Unt = fftn(Un);
   Unft = filter_s .* Unt;
   Unf = ifftn(Unft);
   
   [M, I] = max(reshape(abs(Unf), 1,n^3));
   marble_coords(j, :) = [X(I), Y(I), Z(I)];
   
   isosurface(X,Y,Z, abs(Unf), 0.4)
   axis([-20 20 -20 20 -20 20]), grid on, drawnow
   title("Filtered Signal")
   pause(0.5);
end


%%

figure (100)
plot3(marble_coords(:, 1), marble_coords(:,2), marble_coords(:, 3)), grid on;
xlabel("X coordinate");
ylabel("Y coordinate");
zlabel("Z coordinate");
title("Marble Position (aka you're welcome fluffy)");

%%
final_marble_coordinate_xyz = marble_coords(end, :)
