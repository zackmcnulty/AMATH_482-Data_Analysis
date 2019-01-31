% AMATH 582/482 HW1

clear all; close all; clc;
load Testdata
%%
L=15; % spatial domain 
n=64; % Fourier modes


x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

total_t = zeros(n,n,n);

% for j=1:20 % for each realization of data
%    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
%    figure(1)
%    close all, isosurface(X,Y,Z,abs(Un), 0.4)
%    axis([-20 20 -20 20 -20 20]), grid on, drawnow
%    title('Unfiltered Signal in Spatial Domain')
%    xlabel('X')
%    ylabel('Y')
%    zlabel('Z')
%    set(gca, 'fontsize', 20)
%    pause(2)
% end

%% Averaging of the spectrum

% we have 20 samples of data: time data
% use these different realizations of the object/data to
% average out zero mean frequency noise.

% sum up all frequency data
for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   total_t = total_t + fftn(Un);
end

% normalize the data and shift in the frequency domain so
% it is centered around zero.

ave = abs(fftshift(total_t)) ./ 20;

% find max value so we can normalize
flat = reshape(ave, 1, n^3);
M = max(abs(flat));

ave = ave / M;

% since we normalized ave, the max must be one.
% find frequency index where max occurs
ind2 = find(ave == 1);

% central frequencies in x,y, and z directions respectively.
xc = Kx(ind2);
yc = Ky(ind2);
zc = Kz(ind2);


%% Filtering of the spectrum
close all; clc;

tau = 0.35; % filtering bandwidth

% create 3D gaussian filter. If you are far from central x,y, OR z
% frequency, your signal will be filtered out. This is in shifted
% frequency space however.
filter = exp(-1*tau * ((Kx - xc).^2 +  (Ky - yc).^2 + (Kz - zc).^2));

% unshift filter as it is currently centered at zero.
filter_s = ifftshift(filter);


% apply this filter to frequency domain at each timestep
% to find marbles location at each step
marble_coords = zeros(20, 3);

figure(2)

for j=1:20 % for each realization of data
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   Unt = fftn(Un);
   Unft = filter_s .* Unt;
   Unf = ifftn(Unft);
   
   [M, I] = max(reshape(abs(Unf), 1,n^3));
   marble_coords(j, :) = [X(I), Y(I), Z(I)];
   
   isosurface(X,Y,Z, abs(Unf), 0.4)
   axis([-20 20 -20 20 -20 20]), grid on, drawnow
   hold on;
   title('Filtered Signal in Spatial Domain', 'fontsize' , 20)
   xlabel('X', 'fontsize' , 20);
    ylabel('Y', 'fontsize' , 20);
    zlabel('Z', 'fontsize' , 20);
    xlim([-12, 12]);
    ylim([-10, 10]);
    zlim([-10, 12])
    set(gca, 'fontsize', 25);
    title('Filtered Signal in Spatial Domain', 'fontsize' , 30)
    set(gcf, 'position', [100, 100, 600, 500]);
    saveas(gcf, strcat('images/marble_iso', num2str(j), '.jpg'));
    
    if j == 20
       figure(4)
       view(0,90) % XY
       saveas(gcf, strcat('images/marble_XY.jpg'));
       view(0,0) % XZ
       saveas(gcf, strcat('images/marble_XZ.jpg'));
       view(90,0) % YZ
       saveas(gcf, strcat('images/marble_YZ.jpg'));
    end
end


%%
final_marble_coordinate_xyz = marble_coords(end, :)

figure (2)
plot3(marble_coords(:, 1), marble_coords(:,2), marble_coords(:, 3), 'ro'), grid on;
hold on;
plot3(marble_coords(:, 1), marble_coords(:,2), marble_coords(:, 3), 'b')
hold on;
txt = strcat('\leftarrow time 20: ', mat2str(final_marble_coordinate_xyz));
text(final_marble_coordinate_xyz(1), final_marble_coordinate_xyz(2), final_marble_coordinate_xyz(3), txt, 'fontsize', 20);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Marble Position in Spatial Domain');
set(gca, 'fontsize', 25);
set(gcf, 'position', [100, 100, 600, 500]);
saveas(gcf, strcat('images/marble_position.jpg'));


