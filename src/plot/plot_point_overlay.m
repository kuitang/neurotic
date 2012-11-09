function [ g ] = plot_point_overlay( img, gmm, PHI )

    % Construct an image map. Assume PHI is ordered by x, y.
    [N, D] = size(PHI);
    [X, Y] = size(img);
    c = jet(gmm.K);
    
    overlay = zeros(X, Y, 3);        

    nzidxs   = gmm.n > 0;
    assert(sum(nzidxs) == gmm.K, 'something is horribly wrong');
    
    k_idx    = find(nzidxs);
    k_invidx = zeros(size(gmm.n));    
    % relabel each true class label to a compact class label
    k_invidx(nzidxs) = 1:gmm.K;
    
    for n = 1:N
        x = PHI(n,1);
        y = PHI(n,2);
        overlay(x,y,:) = c(k_invidx(gmm.s_z(n)),:);
    end
    
    %h = imshow(img);
%    hold on
    g = image(overlay);
    hold on
            
    for kk = 2:gmm.K
        k = k_idx(kk);
        text(gmm.mean(k,1), gmm.mean(k,2), num2str(k), ...
             'FontSize', 30, 'FontWeight', 'bold');
        ee = error_ellipse(gmm.pred_cov(1:2,1:2,k), gmm.mean(k,1:2), 'conf', 0.95);
        set(ee, 'color', 'w', 'linewidth', 2);
    end
    hold off
    
    % Plot the 95% (2*sigma) confidence ellipse    
%     for k = 1:mm.K
%         [eig_vec eig_val] = eig(mm.cov(1:2,1:2,k));
%     end
%     
    %set(g, 'AlphaData', 1);    
%    hold off    

end

