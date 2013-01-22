function [ h ] = plot_graph_edge_dist( h, k, mdp, x_max, y_max )
% Visualize the fifth feature, i.e. graph edge distance (based on
% Dijkstra and superpixels)

    % This will show results from the most recent k
    figure(h);
        
    subplot(1, mdp.n_clusters, k);

    dist = mdp.X(:,5);
    dist_img = reshape(dist, x_max, y_max);

    % Transform the column major to row major
    dist_img = dist_img';
    imshow(dist_img, [min(dist) max(dist)]);

    hold on;            
    scatter(mdp.cluster_likes{k}.pdfs{1}.post_mean(1), ...
            mdp.cluster_likes{k}.pdfs{1}.post_mean(2), 'filled');
    hold off;                   
end
