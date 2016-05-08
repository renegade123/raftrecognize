function Coeff = SimulOMP(S, Phi, sigma, T, normType)
% Simultaneous OMP, specify the type of norm used for selecting the atoms
% based on paper (normType = 1 in Tropp's algorithm)
% "Algorithms for Simultaneous Sparse Approximation Part I: Greedy  Pursuit"
% J. Tropp, A. Gilbert, and M. Strauss
% min ||Coeff||_{row,0}   subject to: S = Phi * Coeff
% Input:       S -- signal to be approximated, d x K matrix
%            Phi -- dictionary
%          sigma -- error tolerance
%              T -- number of iterations
%       normType -- type of norm, 1 for L1, 2 for L2, 3 for infinity norm
% Ouptut:  Coeff -- coefficients, a sparse matrix with only T nonzero rows
%         indSet -- index set for selected atoms (common support)
%         Approx -- approximation of S
%            Res -- residuals

% normalizing columns of the dictionary
%Phi0 = Phi;
%Phi = zeros(size(Phi0));
N = size(Phi,2); % the number of atoms
d = size(S,1); K = size(S,2);
Coeff = zeros(N,K); 
Approx = zeros(d,K);

Res = S;
indSet = zeros(T,1);

iter = 1;
norm_res = zeros(T,1);
while ((norm(Res(:)) > sigma) && (iter<=T))
    % compute the projection
    if normType == 1
        [val idx] = sort(sum(abs(Phi'*Res),2), 'descend');
    elseif normType == 2
        [val idx] = sort(sum(abs(Phi'*Res).^2,2), 'descend');
    else
        [val idx] = sort(max(abs(Phi'*Res), [], 2), 'descend');
    end
    
    % update the index set
    indSet(iter) = idx(1);
    
    % 
    Coeff(indSet(1:iter),:) = pinv(Phi(:,indSet(1:iter)))*S;
    
    Approx = Phi*Coeff;
    Res = S - Approx;
    
    norm_res(iter) = norm(Res(:));
    iter = iter + 1;
end

indSet = indSet(1:iter-1);