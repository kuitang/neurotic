function [ mdp ] = refit_graph_dist( mdp, x_max, y_max, k )
% dist = graph_dist(X, G, s) computes Dijkstra vector for shortest paths
% from k's cluster mean of feature matrix X.
    
    if k > 1
        s_x = round(mdp.cluster_likes{k}.pdfs{1}.pred_mean(1));
        s_y = round(mdp.cluster_likes{k}.pdfs{1}.pred_mean(2));

        % Clamp to the image boundary    
        s_x = max(min(1, s_x), x_max);
        s_y = max(min(1, s_y), y_max);

        s = sub2ind([x_max y_max], s_x, s_y);                      

        
        % Compute superpixel distances
        slic_dist = graphshortestpath(mdp.misc_data.edge_G, mdp.misc_data.segments(s), 'Directed', false)';
        % Censor Inf
        slic_dist = min(10000, slic_dist);        
        
        dist = zeros(mdp.N, 1);
        
        % Copy the slic distances to the full distances
        for n_slic = 1:length(slic_dist)
            dist(mdp.misc_data.slic_inds{n_slic}) = slic_dist(n_slic);
        end
        
        mdp.X(:,5) = dist;
        
        % Explicit fit
        kidxs = mdp.cluster_assigns == k;
        mdp.cluster_likes{k}.pdfs{2}.fit(mdp.X(kidxs,5));
        
    end
end
