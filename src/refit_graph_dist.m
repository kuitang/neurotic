function [ mdp ] = refit_graph_dist( mdp, x_max, y_max, k )
% dist = graph_dist(X, G, s) computes Dijkstra vector for shortest paths
% from k's cluster mean of feature matrix X.
    %assert(length(mdp.misc_data.slic_inds{end}) > 1);    
    
    if k > 1
        s_x = round(mdp.cluster_likes{k}.pdfs{1}.post_mean(1));
        s_y = round(mdp.cluster_likes{k}.pdfs{1}.post_mean(2));

        % Clamp to the image boundary
        
        if s_x < 1 || s_x > x_max || s_y < 1 || s_y > y_max
            warning(['crossed image boundary: s_x = ' num2str(s_x) ...
                     ' s_y = ' num2str(s_y)]);
        end
        
        s_x = min(max(1, s_x), x_max);
        s_y = min(max(1, s_y), y_max);                

        s = sub2ind([x_max y_max], s_x, s_y);
        
        % Compute superpixel distances
        my_slic = mdp.misc_data.segments(s);
        [slic_dist, slic_path, ~] = graphshortestpath(mdp.misc_data.edge_G, my_slic, 'Directed', false);
        assert(all(~isinf(slic_dist(:))));
        
        dist = zeros(mdp.N, 1);                
        
        % Copy the normalized slic distances to the full distances
        for n_slic = 1:length(slic_dist)    
            %dist(mdp.misc_data.slic_inds{n_slic}) = slic_dist(n_slic);

            dist(mdp.misc_data.slic_inds{n_slic}) = slic_dist(n_slic) / length(slic_path{n_slic});
        end
        
        % But for all points in our same slic, set distance to zero. Well,
        % not really zero, since we'll evaluate with a Gamma (exponential)
        % distribution.
        dist(mdp.misc_data.slic_inds{my_slic}) = 1e-10;
                                
        % Explicit fit
        kidxs       = mdp.cluster_assigns == k;

        % Shove everything in X, but only fit to the small ones.
        mdp.X(:,5) = dist;        
        mdp.cluster_likes{k}.pdfs{2}.fit(dist(kidxs));
    end
end
