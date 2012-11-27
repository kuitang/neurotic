function [ cluster_assigns ] = initialize_one_cluster_and_background( bg_prob, bg_cutoff, X )
% Produces a cluster initialization for MDP with cluster 1 background,
% cluster 2 everything else. Stupid method sets all points below bg_cutoff
% to background with probability bg_prob.

    [N, D] = size(X);       
    
    % default assign to cluster 2
    cluster_assigns = 2*ones(N, 1);
    
    % Candidate background points
    bgidxs = X(:,3) < bg_cutoff;
    cluster_assigns(bgidxs) = 1;
    
    % And flip everything with probability 1 - bg_prob
    flipidxs = rand(N, 1) < (1 - bg_prob);
    
    cluster_assigns(flipidxs & cluster_assigns == 1) = 2;
    cluster_assigns(flipidxs & cluster_assigns == 2) = 1;

end

