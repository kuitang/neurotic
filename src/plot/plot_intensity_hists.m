function [ h ] = plot_intensity_hists( h, mdp, X )
% Plot intensity histogram for each class.   
    
    K = mdp.n_clusters;
    for k = 1:mdp.n_clusters
        figure(h);
        subplot(2, K, k);
        idxs = mdp.cluster_assigns == k;
        x = X(idxs, 3);
        
        % Top row plots the empirical histograms
        hist(x);
        title(['n = ' num2str(sum(idxs)) ...
               'm = ' num2str(mean(x), 2) ...
               'sd = ' num2str(std(x, 1), 2)]);
           
        % Bottom row plots the predictive distributions
        subplot(2, K, K + k);        
        % DAMN THE SPECIAL CASING!
        if k == 1
            ezplot(@(x) mdp.cluster_likes{1}.pdfs{2}.pred_like_scalar(x), eps, 1);
        else
            m   = mdp.cluster_likes{k}.pred_mean(3);
            sd  = sqrt(mdp.cluster_likes{k}.pred_cov(3, 3));
            dof = mdp.cluster_likes{k}.pred_dof;
            ezplot(@(x) tpdf( (x - m) / sd, dof ), 0, 1);
        end
    end    
end
