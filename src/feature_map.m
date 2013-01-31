function [ PHI, misc_data ] = feature_map( img )
    [X, Y] = size(img);
    N = X * Y;
    
    % Number of features
    D = 5;   
    randon_angles = 180;
    canny_thresh = 0.25;
    
    PHI = zeros(N, D);
    
    bg_intx = radonLikeFeatures(img, true, randon_angles, canny_thresh, @(scanseg) min(scanseg));
    m_bg_intx = mean(bg_intx,3);        
    
    % Scale to 0-1
    % NO, DON'T DO THIS! 
    range = max(m_bg_intx(:)) - min(m_bg_intx(:));
    m_bg_intx = (m_bg_intx - min(m_bg_intx(:))) ./ range;
    %figure; imshow(m_bg_intx);
        
    % Argument 2 is use_raw_imo. This turns on Eq (4) in Kumar.
    bd_intx = radonLikeFeatures(img, false, randon_angles, canny_thresh, @(scanseg) mean(scanseg));
    m_bd_intx = mean(bd_intx,3);
    
    range = max(m_bd_intx(:)) - min(m_bd_intx(:));
    m_bd_intx = (m_bd_intx - min(m_bd_intx(:))) ./ range;
    %figure; imshow(m_bd_intx);
    
    misc_data = struct();    
    misc_data.m_bg_intx = m_bg_intx;
    misc_data.m_bd_intx = m_bd_intx;
    
    misc_data.bd_intx = bd_intx;
    
    misc_data.img = img;
    
    [misc_data.segments, misc_data.slic_inds, misc_data.edge_G] = make_diff_graph(m_bd_intx, 10, 0.1);    
    
    misc_data.slic_dist = graphallshortestpaths(misc_data.edge_G, 'Directed', false);    
    
    n = 1;
    for y = 1:Y
        for x = 1:X            
            PHI(n,:) = [x y img(x,y) m_bg_intx(x,y) double(misc_data.segments(x,y))];
            n = n+1;
        end
    end
    
    % Censor zero intensities (causes some problems for Gamma)
    % Hmmmmmmmm
    PHI(:,3) = max(eps, PHI(:,3));
    PHI(:,4) = max(eps, PHI(:,4));  
end
