function [A_hat E_hat iter] = gzBlockrpcaWithestimatedBlockinfo(D, lambda, tol, maxIter)
% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA.
%
% D - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + lambda * |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end

addpath PROPACK;

% addpath
addpath oriimagedata ;
addpath blockdata ;


%% define block's path
currentPath = cd;
blockPath = fullfile(currentPath,'blockdata') ; % path to files containing initial feature coordinates
userName  = 'smtreeman97320' ;
BlockDir= fullfile(blockPath, userName);
%%

[m n] = size(D);

if nargin < 2
    %lambda = 1 / sqrt(m);
    lambda = 1 / sqrt(m);
end

if nargin < 3
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
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

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;

%gaozhi add
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);
nSetBlock=8;

nWblock=nW/nSetBlock;
nHblock=nH/nSetBlock;

while ~converged       
    iter = iter + 1;
    %% this part for E (sparse matrix) updating; 
    %% i try to enforce the block constraint
    %original code
    %temp_T = D - A_hat + (1/mu)*Y;
    %E_hat = max(temp_T - lambda/mu, 0);
    %E_hat = E_hat+min(temp_T + lambda/mu, 0);
    %%
    %% i modify to enforce the block constraint
    temp_T = D - A_hat + (1/mu)*Y;
    for j = 1:n,
        infName = sprintf('humanBlock%d.mat',j);
        infNames = fullfile(BlockDir, infName);
        load(infNames,'blockInf');
        [nRowb, nColb] = size(blockInf);
        nBlockNum = 0;
        for nB = 1:2:nRowb
            if blockInf(nB+1,2) == 0
               continue;
            else
               nBlockNum=nBlockNum+1;
               Lupx=blockInf(nB,1);
               Lupy=blockInf(nB,2);
               Rdownx=blockInf(nB+1,1);
               Rdowny=blockInf(nB+1,2);
               BlockRect(nBlockNum,1)=Lupx;
               BlockRect(nBlockNum,2)=Lupy;
               BlockRect(nBlockNum,3)=Rdownx;
               BlockRect(nBlockNum,4)=Rdowny;
            end         
        end  
        
        frameTemp = reshape(temp_T(:,j),nH,nW);
        frameCopy = frameTemp;
        for jh = 1:nHblock,
        for jw = 1:nWblock,
            nhStart=(jh-1)*nSetBlock+1;
            nwStart=(jw-1)*nSetBlock+1;
            BlockTemp = frameTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1);
            Blocknorm = norm(BlockTemp,'fro')/sqrt(nSetBlock*nSetBlock);
            
            setratio= 1;
            if Blocknorm < setratio*lambda/mu
               normshrink = 0;
               BlockTemp=BlockTemp*normshrink;
            else 
               normshrink = (Blocknorm - setratio*lambda/mu);
               BlockTemp=BlockTemp*normshrink/Blocknorm;
            end  
            frameTemp(nhStart : nhStart + nSetBlock-1, nwStart : nwStart + nSetBlock-1)=BlockTemp;
        end
        end
        
        [nRowb1, nColb1] = size(BlockRect);
        
        for nB1 = 1:nRowb1
            Lupy=BlockRect(nB1,1)-2;
            Lupx=BlockRect(nB1,2)-2;
            Rdowny=BlockRect(nB1,3)+2;
            Rdownx=BlockRect(nB1,4)+2;
            
            BlockTemp1 = frameCopy(Lupy: Rdowny, Lupx : Rdownx );
            Blocknorm1 = norm(BlockTemp1,'fro')/sqrt((Rdowny - Lupy)*(Rdownx-Lupx));
            
            setratio= 0.000001;
            if Blocknorm1 < setratio*lambda/mu
               normshrink1 = 0;
               BlockTemp1=BlockTemp1*normshrink1;
            else 
               normshrink1 = (Blocknorm1 - setratio*lambda/mu);
               BlockTemp1=BlockTemp1*normshrink1/Blocknorm1;
            end
            
            frameTemp(Lupy: Rdowny, Lupx : Rdownx )=BlockTemp1; 
        end  
        
        frameTemp = reshape(frameTemp,nW*nH,1);
        E_hat(:,j) = frameTemp;
    end
    
    %im1 = uint8(im1);
    %outputFileName = sprintf('A%d.jpg',j);
    %imwrite(im1, outputFileName );   

    if choosvd(n, sv) == 1
        [U S V] = lansvd(D - E_hat + (1/mu)*Y, sv, 'L');
    else
        [U S V] = svd(D - E_hat + (1/mu)*Y, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';    

    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat;
    
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
    end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
