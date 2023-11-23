function [V,explained] = Becket_pca(X,C,lamda)
% [V] = pcaConstrained(X,C)
%   takes an m x n data matrix X and an m x k nuisance matrix C finds
%   the eigenvectors (V) of X, subject to the constraint that these
%   eigenvectors are orthogongonal to the nuisance matrix columns
%
%   because the method uses regression to identify the relationship
%   between X and C, we assumes columns of C are normally distributed
%
% [V,D] = pcaConstrained(X,C) also returns D, a diagonal matrix of
%   eigenvalues corresponding to each eigenvector in V
%
% [V,D,explained] = pcaConstrained(X,C) also returns variance explained
%   by each eigenvector in V; the last element in explained is the variance
%   explained by the nuisance matrix C

    arguments
        X double
        C double
        lamda double
    end

    if size(X,2) ~= size(C,2) % number of columns?
        error("Number of channels won't match")
    end

    % do a little renaming
 

    % extract some variables
    nObs = size(X,1);   
    nObsC = size(C,1);% m observations (X, C rows)
    nFeatures = size(X,2);      % n features (X columns)
    nPCs = min(nObs,nFeatures);
    % preallocate the outputs
    V = NaN(nFeatures,nPCs);    % matrix of eigenvectors

    % center our matrices
    X = X - repmat(nanmean(X,1),nObs,1);
    C = C - repmat(nanmean(C,1),nObsC,1);

    % calculate the covariance of X
    covX = (X'*X) ./ nObs;
    CovC = (C'*C) ./ nObsC;

    % we now need to figure out the relationship between X and C
    % to do this, we borrow from targeted dimensionality reduction
    
    % set up the optimization problem - depend upon type of optimization
    prob = optimproblem('ObjectiveSense','max');
    x = optimvar('x',nFeatures);
    x0.x = randn(nFeatures,1);  % initialize at random
    prob.Objective = (x'*covX*x)/((x'*CovC*x)^lamda); % objective for PCA

    % establish the initial constraints
    prob.Constraints.norm = x'*x == 1; % normal vectors
    
    opts = optimoptions(prob,'Display','final');
    % show(prob) % uncomment to inspect the problem

    % solve for the first PC
    sol = solve(prob,x0,'Options',opts);
    V(:,1) = sol.x;    % we now have the first colum of our eigenvector matrix

    % then iterate through the rest, adding an orthogonality constraint each time
    warning('off','MATLAB:nearlySingularMatrix');
    warning(''); % but first clear last warning so we can respond if the fit goes off the rails
    for k = 2:nPCs
        prob.Constraints.(['ortho',num2str(k-1)]) = x'*V(:,k-1) == 0;
        sol = solve(prob,x0,'Options',opts);
        V(:,k) = sol.x;    % we now have the first colum of our eigenvector matrix
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            warning('Dimensionality of X is < requested nPCs, returning fewer PCs.')
            V = V(:,1:k-1);
            nPCs = k-1;
            break;
        end
    end
    % show(prob) % uncomment to inspect the problem

    % now we also want the eigenvalues to go along with our vectors
    D = V'*covX*V;    % this will give us something reasonable, approach from demixPCA
    D = D.*eye(nPCs); % we'll clean it up though since we have some numerical precision stuff

    % however, to get variance explained, we need to know the total
    % variance, which we won't know from our method b/c we aren't handling
    % nPCs in any kind of detailed way;

    % so we'll approximate our upper bound on variance explained with PCA
    [~,eD] = eig(covX); totalVar = sum(diag(eD));
    % eD = eV'*covX*eV; % note that this gives us eD from eV, c.f. line 75

    % finally, we can return variance explained
    explained = diag(D./totalVar);

    % to figure out the variance explained by nuisance factors
    % we take advantage of the fact that we know sum(explained) == 1
    explained(end+1) = 1-sum(explained);
    % CAUTION: THIS ONLY WORKS IF NUM PCs is == the dimensionality of the
    % matrix after controlling for the nuisance variables

end