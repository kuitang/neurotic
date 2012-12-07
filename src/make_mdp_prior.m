function [ mdp ] = make_mdp_prior( X )
    [N, D] = size(X);
    nX = max(X(:,1));
    nY = max(X(:,2));
    
    %% Parameter priors
    prior_mean = [100 100 0.75 0.4]; % mean of the nonzero entries
    
    prior_cov  = cov(X);
    
    v_intensity = prior_cov(3,3);
    v_bg_radon  = prior_cov(4,4);
    
    prior_cov(3:4,:) = zeros(2, 4);
    prior_cov(:,3:4) = zeros(4, 2);    
        
    prior_cov(3,3) = 0.5 * v_intensity;
    prior_cov(4,4) = 0.5 * v_bg_radon;
    
    background = GammaGamma(1, 6, 2);
    
    prior_dof = 5;
    prior_n   = 0.5;
    
    concentration = 1;        
    
    %% Initialize cluster assignments
    % Wood: assign them uniformly at random
    cluster_assigns = initialize_one_cluster_and_background(0.99, 0.5, X);
    
    %% Construct objects
    nw_prior = NormalWishart(prior_mean, prior_cov, prior_dof, prior_n);
    mdp = MDP(concentration, X, cluster_assigns, nw_prior);
    
    
    % Hack: one uniform for (x,y) and another uniform for dimension 4 (the
    % background Radon response) ... ; NO IT'S BORKED!
    mdp.cluster_likes{1} = ProductDistribution(1, 2, Uniform(size(X, 1)), ...
                                               3, 3, background, ...
                                               4, 4, Uniform(1)); % Because domain is [0,1]
    mdp.refit(1);
    %[particles, weights] = smc_init(mdp, ProductDistribution(1, 2, Uniform(size(X, 1)), 3, 3, background));
    
    % We will take the MODE assignment. (Wait, fuck it, the model need not
    % be identifiable...)

    % Add the background class
    % TODO: Change or sensitivity analysis!            
end
