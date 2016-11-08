function RelativeChange = gzComputeRelativeChange(data, A_hat)
%[A_hat E_hat iter] = inexact_alm_rpca(D, lambda, tol, maxIter)

[m n] = size(data);
RelativeChange = zeros( m, n);

Change         =  abs(data - A_hat);
RelativeChange =  Change./data;

end
