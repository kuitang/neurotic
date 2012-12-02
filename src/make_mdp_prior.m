function [ mdp ] = make_mdp_prior( X )
    [N, D] = size(X);
    nX = max(X(:,1));
    nY = max(X(:,2));
    
    %% Parameter priors
    prior_mean = [100 100 0.75];
    
    prior_cov  = cov(X);
    v_intensity = prior_cov(3,3);
    prior_cov(3,:) = zeros(3, 1);
    prior_cov(:,3) = zeros(1, 3);   
    prior_cov(3,3) = 0.1 * v_intensity;
    
    background = GammaGamma(1, 6, 2);
    
    prior_dof = 5;
    prior_n   = 0.5;
    
    concentration = 1;        
    
    %% Initialize cluster assignments
    % Wood: assign them uniformly at random
    cluster_assigns = initialize_one_cluster_and_background(0.99, 0.4, X);
    
    %% Construct objects
    nw_prior = NormalWishart(prior_mean, prior_cov, prior_dof, prior_n);
    mdp = MDP(concentration, X, cluster_assigns, nw_prior);
    
    mdp.cluster_likes{1} = ProductDistribution(1, 2, Uniform(size(X, 1)), 3, 3, background);
    mdp.refit(1);
    %[particles, weights] = smc_init(mdp, ProductDistribution(1, 2, Uniform(size(X, 1)), 3, 3, background));
    
    % We will take the MODE assignment. (Wait, fuck it, the model need not
    % be identifiable...)

    % Add the background class
    % TODO: Change or sensitivity analysis!            
end
