function [ gmm ] = test_gmm(niter)
% Sanity check your algorithm with iris
    
    X = iris_dataset';    
    [N D] = size(X);
    gmm = struct('K', 3, ...
                 'prior_mean', mean(X), ...
                 'prior_cov', eye(D) / 2, ... % Because expectation of Wishart is dof * c.
                 'prior_scale', 1, ...
                 'prior_mix', 10, ...
                 'prior_dof', 5, ...
                 's_z', randsample(3,N,true));                 
    
    gmm = gmm_gibbs(X, [], gmm, niter, 10, 50);
end
