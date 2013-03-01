function [ mdp, state, params ] = make_mdplight_prior( X, misc_data )
    [N, D] = size(X);
    nX = max(X(:,1));
    nY = max(X(:,2));
    
    %% Parameter priors
    prior_mean = [100 100 0.9 0.4]'; % mean of the nonzero entries
    prior_cov  = cov(X(:,1:4));
    
    v_intensity = prior_cov(3,3);
    v_bg_radon  = prior_cov(4,4);
    
    prior_cov(3:end,:) = 0;
    prior_cov(:,3:end) = 0;    
        
    prior_cov(3,3) = 0.5 * v_intensity;
    prior_cov(4,4) = 0.5 * v_bg_radon;        
    
    prior_dof = 5;
    prior_n   = 2;
    
    concentration = 1;        
        
    prior = struct('rate_shape', 20, 'rate_rate', 10, ...
                   'mean', prior_mean, 'n', prior_n, ...
                   'dof', prior_dof, 'prec_chol', chol(inv(prior_cov)), ...
                   'bg_shape', 10, 'bg_rate', 5);
    
    % cluster_assigns specifies two clusters. Start the Markov chain at
    % their priors.
    %
    % Initialization is really arbitrary...
    state(1).bg_rate = 10;    
    state(2).mean = prior_mean;
    state(2).rate = 10;
    state(2).prec_chol = chol(inv(prior_cov));  
    
    cluster_assigns = initialize_one_cluster_and_background(0.99, 0.7, X);               
    mdp = MDPLight(concentration, X, prior, state, misc_data, cluster_assigns);   
    
    graph_dist = GraphDist(misc_data);
    
    params = struct('m', 3, 'b', 3, ...
                    'precision_slice_sample', false, 'mean_slice_sample', false, ...
                    'bg_Z', length(X), 'graph_dist', graph_dist, ...
                    'mu_prop_chol', chol(diag([1/20 1/20 1/5 1/5])));


    % Add the background class
    % TODO: Change or sensitivity analysis!            
end
