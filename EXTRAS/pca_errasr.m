function [V,explained] = pca_errasr(X,X_stim,k,lamda)
% Sparse principal component analysis based on optimization over Stiefel.
%
% [V,explained] = pca_errasr(X,X_stim,lamda ,NameValueArgs)
%
% Here We will compute principal components based on the cost function 
% U.T cov(X_stim) U / (U.T cov(X) U)**lamda - maximize over stiefel
% manifold ( orthogonal components constraines)
%
% It will p bring the variance of stim more together and push the neural
% variance down - possibly good components to filter out stim artifacts
% 
% [V,explained] = pcaConstrained(X,C) also returns variance explained
%   by each eigenvector in V; the last element in explained is the variance
%   explained by the nuisance matrix C

    if ~exist('lamda', 'var') || isempty(lamda)
        lamda = 0.05;
    elseif ~exist('k', 'var') || isempty(k)
        k = size(X_stim,2);
    elseif nargin <2
        error('Please provide 3 inputs (or none for a demo).');
    end
    
    % do the covariancce and mean centerise the matrix
    X_centered = X - mean(X, 1);
    Xstim_centered = X_stim - mean(X_stim, 1);


    [n, p] = size(X_centered);

    St = stiefelfactory(p, k); %this might need to be configured
    problem.M = St;
  
    % Define the cost function here and set it in the problem structure.
    problem.cost = @cost;

    problem.egrad = @egrad;
    
    [V, ~, info, ~] = trustregions(problem);
    semilogy([info.iter], [info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the gradient of f');

    % Project the data onto the top k principal components
    projected_data = Xstim_centered*V;
    
    % Variance of the projected data
    variance_projected = sum(var(projected_data, 0, 2));
    
    % Total variance of the original data
    cov_matrix = cov(Xstim_centered);
    total_variance = trace(cov_matrix);
    
    % Calculate the fraction of variance explained
    explained = variance_projected / total_variance;
    
    function G = egrad(V)
        t0 = norm(X_centered*V, 'fro');
        t1 = norm(Xstim_centered*V, 'fro');
        a = 2*lamda;
        cov_xstim = (Xstim_centered' * Xstim_centered); % not really a covariance
        cov_x = (X_centered' * X_centered);

        G1 = (2*(t0^(-a))) *(cov_xstim *V);

        G2 = a*(t0^(-(1+a)))*((t1^2)/(t0))*cov_x*V;
        G = -(G1-G2);
    end

    function val = cost(V)
        val = -norm(Xstim_centered*V, 'fro')^2/(norm(X_centered*V, 'fro')^(2*lamda));
    end

    
end
