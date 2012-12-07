function [ h ] = plot_img( h, n, loglike_vec, X, img, mdp )
% Plot images for an mdp.
% h           - handle to plot into
% n           - iteration number
% loglike_vec - vector of log likelihoods. Trailing zeros ignored.
% img         - matrix of the image
% mdp         - MDP, or really anything with .cluster_assigns
    
    figure(h);
    
    subplot(1,4,1);
    nz_loglike = loglike_vec(loglike_vec ~= 0);
    plot(nz_loglike);
    title('Log-likelihood vs iteration');
    xlabel('Iteration');
    ylabel('Log-likelihood');

    subplot(1,4,2);
    plot_point_overlay(img, mdp, X);
    title(['K = ' num2str(mdp.n_clusters) ' iter ' num2str(n) ' ll = ' num2str(nz_loglike(end))]);

    subplot(1,4,3);
    imshow(img);
    
    % The feature map
    subplot(1,4,4);
    fmap = reshape(X(:,4), size(img));
    imshow(fmap);

end
