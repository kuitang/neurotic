function [ mdp ] = make_mdp_prior( X, misc_data )
    [N, D] = size(X);
    nX = max(X(:,1));
    nY = max(X(:,2));
    
    %% Parameter priors
    prior_mean = [100 100 0.9 0.4]; % mean of the nonzero entries
    
    prior_cov  = cov(X(:,1:4));
    
    v_intensity = prior_cov(3,3);
    v_bg_radon  = prior_cov(4,4);
    
    prior_cov(3:end,:) = 0;
    prior_cov(:,3:end) = 0;    
        
    prior_cov(3,3) = 0.5 * v_intensity;
    prior_cov(4,4) = 0.5 * v_bg_radon;    
    
    background = GammaGamma(1, 10, 5);
    
    prior_dof = 5;
    prior_n   = 0.5;
    
    concentration = 1;        
    
    %% Initialize cluster assignments
    % Wood: assign them uniformly at random
    cluster_assigns = initialize_one_cluster_and_background(0.99, 0.7, X);
    
    %% Construct objects
        
    feature_prior  = NormalWishart(prior_mean, prior_cov, prior_dof, prior_n);
    distance_underlying = GammaGamma(1, 20, 10);
    
    
    % HACK; to patch up later.
    class_prior = ProductDistribution(1, 4, feature_prior, ...
                                      5, 5, distance_underlying);
    
    % Hack: low probabilities when plugged into NormalGamma.
    % Interpretation: when starting a new cluster, take its Dijkstra
    % distance to be "far".
                                          
    mdp = MDP(concentration, X, misc_data, cluster_assigns, class_prior);    
    
    distance_prior = GraphDistDistribution(distance_underlying);
    distance_post  = GraphDistDistribution(mdp.cluster_likes{2}.pdfs{2});
        
    distance_prior.mdp = mdp;
    distance_post.mdp = mdp;
    distance_prior.misc_data = misc_data;        
    distance_post.misc_data = misc_data;
    
    distance_post.my_k = 2;
    
    % Hack: plug back in the distance underlying. mdp.prior is a
    % ProductDistribution
    mdp.prior.pdfs{2} = distance_prior;
    mdp.cluster_likes{2}.pdfs{2} = distance_post;        
    
    % Hack: one uniform for (x,y) and another uniform for dimension 4 (the
    % background Radon response)
    mdp.cluster_likes{1} = ProductDistribution(1, 2, Uniform(size(X, 1)), ...
                                               3, 3, background, ...
                                               4, 4, Uniform(1), ... % Domain [0, 1]
                                               5, 5, Uniform(1)); % Domain essentially 1        
    
    %[particles, weights] = smc_init(mdp, ProductDistribution(1, 2, Uniform(size(X, 1)), 3, 3, background));
    
    % We will take the MODE assignment. (Wait, fuck it, the model need not
    % be identifiable...)

    % Add the background class
    % TODO: Change or sensitivity analysis!            
end
