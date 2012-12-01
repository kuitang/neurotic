function [ PHI ] = feature_map( img )
    [X, Y] = size(img);
    N = X * Y;
    
    D = 3;
    
    PHI = zeros(N, D);
    n = 1;
    for y = 1:Y
        for x = 1:X
            PHI(n, :) = [x y img(x,y)];
            n = n+1;
        end
    end
    
    % Censor zero intensities (causes some problems for Gamma)
    PHI(:,3) = max(eps, PHI(:,3));
end

