function [A_hat,E_hat,dt_hat,numIter] = rasl_inner_apg(D, J, lambda,  ...
    tol, maxIter, continuationFlag, mu )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input :
%     D - m x n matrix of observations/data (required input)
%     J - m x pt matrix of the Jacobian matrix, pt is the number of
%         parameters of the tranformation
%     lambda - weight on sparse error term in the cost function (required input)
%     
%     tol - tolerance for stopping criterion.
%         - DEFAULT 1e-6 if omitted or -1.
%     maxIter - maximum number of iterations
%             - DEFAULT 10000, if omitted or -1.
%     continuationFlag - 1 if a continuation is to be done on the parameter mu
%                      - DEFAULT 1, if omitted or -1.
%     mu - relaxation parameter
%        - ignored if continuationFlag is 1.
%        - DEFAULT 1e-3, if omitted or -1.
%     outputFileName - Details of each iteration are dumped here, if provided.

% output :
%     [A_hat, E_hat, dt_hat] - estimates for the low-rank part, error part and transformation, respectively
%     numIter - number of iterations until convergence
%     numIter - the number of iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
    error('Too few arguments') ;
end

if nargin < 4
    tol = 1e-6 ;
elseif tol == -1
    tol = 1e-6 ;
end

if nargin < 5
    maxIter = 1000 ;
elseif maxIter == -1
    maxIter = 1000 ;
end

if nargin < 6
    continuationFlag = 1 ;
elseif continuationFlag == -1 ;
    continuationFlag = 1 ;
end

if nargin < 7
    mu = 1e-3 ;
elseif mu == -1
    mu = 1e-3 ;
end

DISPLAY_EVERY = 100 ;

%% Initializing optimization variables

[m,n] = size(D) ;

maxJNorm = norm(J{1}) ;
for fileIndex = 2 : n
    maxJNorm = max(maxJNorm,norm(J{fileIndex})) ;
end
tau_0 = sqrt(3*(2+maxJNorm*maxJNorm)*max(1,maxJNorm*maxJNorm)) ;  
% Lipschitz constant

X_km1_A = zeros(m,n) ; X_km1_E = zeros(m,n) ; X_km1_dt = cell(1,n) ; % X^{k-1} = (A^{k-1},E^{k-1}, dt^{k-1})
X_k_A = zeros(m,n) ; X_k_E = zeros(m,n) ; X_k_dt = cell(1,n) ; % X^{k} = (A^{k},E^{k},E^{k})
Jn = size(J{1}, 2) ; % the number of parameters of the transformation
for i = 1:n
    X_km1_dt{i} = zeros(Jn, 1) ; 
    X_k_dt{i} = zeros(Jn, 1) ;
end

Delta = zeros(m,n) ; % auxilliary variable 
delta_xi = cell(1,n) ;

if continuationFlag
    mu_0 = norm(D) ;
    mu_k = mu_0 ;
    mu_bar = 1e-4 * mu_0 ;    
else
    mu_k = mu ;
end

tau_k = tau_0 ;

t_k = 1 ; % t^k
t_km1 = 1 ; % t^{k-1}

converged = 0 ;
numIter = 0 ;

%% Start main inner loop

while ~converged
    
    tt = (t_km1 - 1)/t_k ;
    Y_k_A = X_k_A + tt*(X_k_A-X_km1_A) ;
    Y_k_E = X_k_E + tt*(X_k_E-X_km1_E) ;
    for i = 1 : n
        Y_k_dt{i} = X_k_dt{i} + tt*(X_k_dt{i}-X_km1_dt{i}) ;
    end
    
    for i = 1 : n
        Delta(:,i) = J{i}*Y_k_dt{i} ;
    end
    G_k = Y_k_A+Y_k_E-D-Delta ;
    G_k_A = Y_k_A - (1/tau_k)*G_k ;
    G_k_E = Y_k_E - (1/tau_k)*G_k ;
    for i = 1 : n
        delta_xi{i} = - J{i}'*G_k(:,i) ;
        G_k_dt{i} = Y_k_dt{i} - (1/tau_k)*delta_xi{i} ;
    end
        
    [U,S,V] = svd(G_k_A,'econ') ;
    diagS = diag(S) ;
        
    X_kp1_A = U * diag(pos(diagS-mu_k/tau_k)) * V';
    X_kp1_E = sign(G_k_E) .* pos( abs(G_k_E) - lambda*mu_k/tau_k );
    X_kp1_dt = G_k_dt ;
        
    rankA  = sum(diagS>mu_k/tau_k);
    cardE = sum(sum(double(abs(X_kp1_E)>0)));
    
    numIter = numIter + 1 ;
    
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    % stopping condition
    for i = 1 : n
        Delta(:,i) = J{i}*(X_kp1_dt{i} - Y_k_dt{i}) ;
    end
    temp = X_kp1_A + X_kp1_E - Y_k_A - Y_k_E - Delta ;
    S_kp1_A = tau_k*(Y_k_A-X_kp1_A) + temp ;
    S_kp1_E = tau_k*(Y_k_E-X_kp1_E) + temp ;
    for i = 1 : n
        delta_xi{i} = - J{i}'*temp(:,i) ;
        S_kp1_dt{i} = tau_k*(Y_k_dt{i}-X_kp1_dt{i}) + delta_xi{i} ;
    end
    stoppingCriterion = (norm([S_kp1_A,S_kp1_E],'fro')+norm(cell2mat(S_kp1_dt),'fro'))/(tau_k*max(1,(norm([X_kp1_A,X_kp1_E],'fro')+norm(cell2mat(X_kp1_dt),'fro')))) ;
    
    if mod(numIter,DISPLAY_EVERY) == 0
        disp(['Iteration ' num2str(numIter) '  rank(A) ' num2str(rankA) ...
            ' ||E||_0 ' num2str(cardE)  '  Stopping Criterion ' ...
            num2str(stoppingCriterion)]) ;
    end
    
    if stoppingCriterion <= tol
        disp(['RASL inner loop is converged at:']);
        disp(['Iteration ' num2str(numIter) '  rank(A) ' num2str(rankA) ...
            ' ||E||_0 ' num2str(cardE)  '  Stopping Criterion ' ...
            num2str(stoppingCriterion)]) ;
        converged = 1 ;
    end
    
    % continuation
    if continuationFlag
        mu_k = max(0.9*mu_k,mu_bar) ;
    end
    
    % update the data
    t_km1 = t_k ;
    t_k = t_kp1 ;
    X_km1_A = X_k_A ; X_km1_E = X_k_E ; X_km1_dt = X_k_dt ;
    X_k_A = X_kp1_A ; X_k_E = X_kp1_E ; X_k_dt = X_kp1_dt ;
        
    if ~converged && numIter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
    
end

% output
A_hat = X_k_A ;
E_hat = X_k_E ;
dt_hat = X_k_dt ;