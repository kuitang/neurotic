function [ h ] = plot_intensity_hists( h, mdp, X )
% Plot intensity histogram for each class.   
    
    K = mdp.n_clusters;
    xs = linspace(eps, 1, 100);
    
    for k = 1:K
        figure(h);
        idxs = mdp.cluster_assigns == k;
        x = X(idxs, 3);
        
        % Top row plots the empirical histograms
        subplot(3, K, k), hist(x);
        title(['n = ' num2str(sum(idxs)) ...
               ' m = ' num2str(mean(x), 2) ...
               ' sd = ' num2str(std(x, 1), 2)]);
                   
        subplot(3, K, K + k);
        
        % DAMN THE SPECIAL CASING!
        if k == 1
            plot(xs, mdp.cluster_likes{1}.pdfs{2}.pred_like_scalar(xs));            
        else
            m   = mdp.cluster_likes{k}.pdfs{1}.post_mean(3);            
            sd  = sqrt(mdp.cluster_likes{k}.pdfs{1}.pred_cov(3, 3));
            dof = mdp.cluster_likes{k}.pdfs{1}.pred_dof;
            plot(xs, tpdf( (xs - m) / sd, dof ));            
        end
        
        % This is plotting the GammaGamma
        subplot(3, K, 2*K + k);
        
        if k ~= 1
            %xxs = linspace(eps, 500, 100);
            %plot(xxs, mdp.cluster_likes{k}.pdfs{2}.underlying.pred_like_scalar(xxs));
            hist(mdp.cluster_likes{k}.pdfs{2}.empirical_distribution());
        end            
        
    end    
end
