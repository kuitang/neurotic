function [ g ] = plot_point_overlay( img, gmm, PHI )

    % Construct an image map
    [N, D] = size(PHI);
    [X, Y] = size(img);
    c = jet(gmm.K);
    
    overlay = zeros(X, Y, 3);        
    
    for n = 1:N
        x = PHI(n,1);
        y = PHI(n,2);
        overlay(x,y,:) = c(gmm.s_z(n),:);
    end
    
    h = imshow(img);
    hold on
    g = image(overlay);    
    set(g, 'AlphaData', 1);    
    hold off    

end

