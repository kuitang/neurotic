function [ g ] = plot_point_overlay( img, mdp, X )
% Plot point assignments and 95% confidence ellipses.
% img - a grayscale image (array) of the original data
% mdp - an MDP whose cluster_likes have a .post_mean and .post_cov
% X   - the data matrix whose first 3 columns are [x y intensity]

    % Construct an image map. Assume PHI is ordered by x, y.
    [N, D] = size(X);
    [nX, nY] = size(img);        
        
    c = jet(mdp.n_clusters);
    
    overlay = zeros(nX, nY, 3);            
    
    for n = 1:N
        x = X(n,1);
        y = X(n,2);
        overlay(x,y,:) = c(mdp.cluster_assigns(n), :);        
    end
    
    g = image(overlay);
    hold on
    
    for k = 2:mdp.n_clusters
        % We'll cheat here and plot the error ellipses of the PREDICTIVE
        % distributions. This is visually accurate: it shows how likely a
        % point is to be classified in a given cluster.
        %
        % A more correct plot, encapsulating the uncertainty of the model
        % itself, would be ellipses of the posterior mean.
        
        cl = mdp.cluster_likes{k};        
        mean_x = cl.pred_mean(1);
        mean_y = cl.pred_mean(2);
                                
        text(mean_x, mean_y, num2str(k), ...
             'FontSize', 30, 'FontWeight', 'bold');
        
        ee = error_ellipse(cl.pred_cov(1:2,1:2), cl.pred_mean(1:2), 'conf', 0.95);
        set(ee, 'color', 'w', 'linewidth', 2);
    end    
    hold off  

end

