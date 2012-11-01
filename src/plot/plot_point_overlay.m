function [ g ] = plot_point_overlay( img, gmm, PHI )

    % Construct an image map. Assume PHI is ordered by x, y.
    [N, D] = size(PHI);
    [X, Y] = size(img);
    c = jet(gmm.K);
    
    overlay = zeros(X, Y, 3);        
    
    for n = 1:N
        x = PHI(n,1);
        y = PHI(n,2);
        overlay(x,y,:) = c(gmm.s_z(n),:);
    end
    
    %h = imshow(img);
%    hold on
    g = image(overlay);
    hold on
    for k = 1:gmm.K
        text(gmm.mean(k,1), gmm.mean(k,2), num2str(k), ...
             'FontSize', 30, 'FontWeight', 'bold');
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
