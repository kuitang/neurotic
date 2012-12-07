function [ PHI ] = feature_map( img )
    [X, Y] = size(img);
    N = X * Y;
    
    % Number of features
    D = 4;
    
    res = 72;
    
    PHI = zeros(N, D);
    
    bg_intx = radonLikeFeatures(img, true, res, @(scanseg) min(scanseg));
    m_bg_intx = mean(bg_intx,3);
    
    % Scale to 0-1
    range = max(m_bg_intx(:)) - min(m_bg_intx(:));
    m_bg_intx = (m_bg_intx - min(m_bg_intx(:))) ./ range;
    
    n = 1;
    for y = 1:Y
        for x = 1:X
            PHI(n, :) = [x y img(x,y) m_bg_intx(x,y)];
            n = n+1;
        end
    end
    
    % Censor zero intensities (causes some problems for Gamma)
    PHI(:,3) = max(eps, PHI(:,3));
    PHI(:,4) = max(eps, PHI(:,4));
  
end
