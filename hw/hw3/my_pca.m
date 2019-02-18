function A = my_pca(rank_approx, varargin)

% rank_approx: What dimension of rank-approximation we want to perform.
% i.e. how many of the singular values are relevant?

x1 = varargin{1};
y1 = varargin{2};
x2 = varargin{3};
y2 = varargin{4};
x3 = varargin{5};
y3 = varargin{6};

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
vars = var(X.').';

% demean data
X = X - means;
%X = (X - means) ./ vars;

[u, s, v] = svd(X);

% By diagonalizing our covariance matrix, we can generate some
% completely independent component (see pg 117 of AMATH 582/482 notes)

Y = u'*X; % note we want ' not .' as here we want complex conjugate U*
principal_components = Y(1:rank_approx, :)

% plot the normalized sigma values
singular_values = diag(s);
plot(singular_values ./ max(singular_values), 'ro');
title("Normalized Singular Values")
ylabel('\sigma_j')
xlabel('index j')


% Reconstructing X using a low-rank approximation
A = zeros(size(X));
whos

for j = 1:rank_approx
    A = A + singular_values(j).*u(:, j)*v(:, j)'; % transpose is our *
end

A = A + means;
%A = A .* vars  + means;
% energy captured?

end

