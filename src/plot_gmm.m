function [ h ] = plot_gmm( img, gmm )
    % h = plot_gmm(img, gmm) overlays covariance circles on img        
    
    h = imshow(img);
    hold on
    c = jet(gmm.K);
    
    min_sz = min(1 ./ gmm.lam(:,3));
    scale = 100 / min_sz;
    
    for k = 1:gmm.K
        x = gmm.mu(k,1);
        y = gmm.mu(k,2);
        s = scale / gmm.lam(k,3);
        scatter(x, y, s, 'r');
    end

end

