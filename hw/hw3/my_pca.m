function A = my_pca(rank_approx, pcs, offset, yrange, varargin)

% rank_approx: What dimension of rank-approximation we want to perform.
% i.e. how many of the singular values are relevant?
% pcs = number of principal components to plot
% yrange = for plotting; y limits on low rank approximations
% varargin = time series data; list of vectors x1, y1, x2, etc...

x1 = varargin{1};
y1 = varargin{2};
x1 = x1(offset(1):end);
y1 = y1(offset(1):end);

x2 = varargin{3};
y2 = varargin{4};
x2 = x2(offset(2):end);
y2 = y2(offset(2):end);

x3 = varargin{5};
y3 = varargin{6};
x3 = x3(offset(3):end);
y3 = y3(offset(3):end);

% add all time measurments to a matrix with time varying across columns 
% and position measurements as rows. i.e.
% [x1(1) x2(2) ....;
%   y1(1) y2(2) ....]
% This makes U in the SVD describe the principal directions in space and
% V the principal directions in time??

% position vectors are not all the same length; find the min length and 
% use that many points instead.
n = min([length(x1); length(x2); length(x3)]);
X = [x1(1:n); y1(1:n); x2(1:n) ; y2(1:n); x3(1:n); y3(1:n)];
means = mean(X.').';


% demean data
X = X - means;


[u, s, v] = svd(X);

% By diagonalizing our covariance matrix, we can generate some
% completely independent component (see pg 117 of AMATH 582/482 notes)

Y = u'*X; % note we want ' not .' as here we want complex conjugate U*
principal_components = Y(1:pcs, :);

% plot the normalized sigma values
singular_values = diag(s);
plot(singular_values ./ max(singular_values), 'r.', 'markersize', 40);
xticks(1:6)
title("Normalized Singular Values")
ylabel('\sigma_j')
xlabel('index j')
set(gca, 'fontsize', 20);



% Reconstructing X using a low-rank approximation
A = zeros(size(X));

for j = 1:rank_approx
    A = A + singular_values(j).*u(:, j)*v(:, j)'; % transpose is our *
end

% Remean the data!
A = A + means;
     

% energy captured?
energy = sum(singular_values(1:rank_approx)) / sum(singular_values)


figure(7)
subplot(231)
plot(x1, 'r'), hold on;
ylim(yrange)
xlim([0, n])
title('x1')
plot(A(1,:), 'b')
legend({'original', 'low-rank approx'})
set(gca, 'fontsize', 15);
ylabel('Position')


subplot(234)
plot(y1, 'r'), hold on
ylim(yrange)
xlim([0, n])
title('y1')
plot(A(2,:), 'b')
legend({'original', 'low-rank approx'})
set(gca, 'fontsize', 15);
ylabel('Position')


subplot(232)
plot(x2, 'r'), hold on;
ylim(yrange)
xlim([0, n])
title('x2')
plot(A(3,:), 'b')
legend({'original', 'low-rank approx'})
set(gca, 'fontsize', 15);


subplot(235)
plot(y2, 'r'), hold on;
ylim(yrange)
xlim([0, n])
title('y2')
plot(A(4,:), 'b')
legend({'original', 'low-rank approx'})
xlabel('Time (frame number)')
set(gca, 'fontsize', 15);



subplot(233)
plot(x3, 'r'), hold on;
ylim(yrange)
xlim([0, n])
title('x3')
plot(A(5,:), 'b')
legend({'original', 'low-rank approx'})
set(gca, 'fontsize', 15);


subplot(236)
plot(y3, 'r'), hold on;
ylim(yrange)
xlim([0, n])
title('y3')
plot(A(6,:), 'b')
legend({'original', 'low-rank approx'})
set(gca, 'fontsize', 15);



% plot the principal_components
figure(8)
names = [];
for k = 1:pcs
   plot(principal_components(k, :), 'linewidth', 2), hold on; 
   names = [names strcat("Principal Component ", num2str(k))];
   xlabel('Time (frame)')
   ylabel('Position')
   
end
   xlim([0, length(principal_components)]);
   title('Principal Components');
   legend(names);
   set(gca, 'fontsize', 20);

end

