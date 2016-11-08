% test lan
% 
%
addpath PROPACK;

% initialize
D=rand(10,3);
Y = D;
norm_two = lansvd(Y, 1, 'L');

norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);

Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
mu = 1.25/norm_two % this one can be tuned
mu_bar = mu * 1e7
rho = 1.5          % this one can be tuned
d_norm = norm(D, 'fro');
