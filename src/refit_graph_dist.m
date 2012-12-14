function [ mdp ] = refit_graph_dist( mdp, x_max, y_max, k )
% dist = graph_dist(X, G, s) computes Dijkstra vector for shortest paths
% from k's cluster mean of feature matrix X.
    %assert(length(mdp.misc_data.slic_inds{end}) > 1);

    if k > 1
        s_x = round(mdp.cluster_likes{k}.pdfs{1}.pred_mean(1));
        s_y = round(mdp.cluster_likes{k}.pdfs{1}.pred_mean(2));

        % Clamp to the image boundary    
        s_x = max(min(1, s_x), x_max);
        s_y = max(min(1, s_y), y_max);

        s = sub2ind([x_max y_max], s_x, s_y);                      

        
        % Compute superpixel distances
        my_slic = mdp.misc_data.segments(s);
        slic_dist = graphshortestpath(mdp.misc_data.edge_G, my_slic, 'Directed', false)';
        % Censor Inf
        cap = 10000;
        slic_dist = min(cap, slic_dist);        
        
        dist = zeros(mdp.N, 1);
        
        %assert(length(mdp.misc_data.slic_inds{end}) > 1);
        
        % Copy the slic distances to the full distances
        for n_slic = 1:length(slic_dist)     
            %assert(length(mdp.misc_data.slic_inds{end}) > 1);            
            dist(mdp.misc_data.slic_inds{n_slic}) = slic_dist(n_slic);        
        end
        
        % But for all points in our same slic, set distance to zero.
        dist(mdp.misc_data.slic_inds{my_slic}) = 0;
                
        %mdp.X(:,5) = dist;
        
        % Explicit fit
        kidxs       = mdp.cluster_assigns == k;
        finite_idxs = dist < cap;

        % Shove everything in X, but only fit to the small ones.
        mdp.X(:,5) = dist;
        works = dist(kidxs & finite_idxs);
                
        mdp.cluster_likes{k}.pdfs{2}.fit(works);
    end
end
