function [ cluster_assigns ] = initialize_clusters_unif_with_background( background_like, X )
% Produces a cluster initialization for MDP.
    [N, D] = size(X);       
    
    % Sample with replacement
    cluster_assigns = randsample(2:N, N, true)'; % Transpose for column
    
    assert(all(cluster_assigns > 0));
    
    % Compactify
    k = 1;
    n = N;
    while k < n
        if ~any(cluster_assigns == k)
            gtidxs = cluster_assigns > k;
            cluster_assigns(gtidxs) = cluster_assigns(gtidxs) - 1;
            n = n - 1;
        else
            k = k + 1;
        end
    end

    bg = rand < background_like(X(:,3));
    % assign to background
    cluster_assigns(bg) = 1;            

end

