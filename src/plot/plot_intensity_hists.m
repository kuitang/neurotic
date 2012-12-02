function [ h ] = plot_intensity_hists( h, mdp, X )
% Plot intensity histogram for each class.   
    
    K = mdp.n_clusters;
    
    for k = 1:K
        figure(h);
        idxs = mdp.cluster_assigns == k;
        x = X(idxs, 3);
        
        % Top row plots the empirical histograms
        subplot(2, K, k), hist(x);
        title(['n = ' num2str(sum(idxs)) ...
               ' m = ' num2str(mean(x), 2) ...
               ' sd = ' num2str(std(x, 1), 2)]);
                   
        subplot(2, K, K + k);
        xs = linspace(eps, 1, 100);
        % DAMN THE SPECIAL CASING!
        if k == 1
            plot(xs, mdp.cluster_likes{1}.pdfs{2}.pred_like_scalar(xs));            
        else
            m   = mdp.cluster_likes{k}.post_mean(3);
            sd  = sqrt(mdp.cluster_likes{k}.pred_cov(3, 3));
            dof = mdp.cluster_likes{k}.pred_dof;
            plot(xs, tpdf( (xs - m) / sd, dof ));            
        end
    end    
end
