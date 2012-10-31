function [ like ] = background_filter(gmm, X, radius)
    like = zeros(gmm.N, 1);
    
    sz   = (1 + 2 * radius)^2;
    mass = sz/2 * (sz + 1);
    
    % TODO: Vectorize
    
    for n = 1:gmm.N
        x = X(n,1);
        y = X(n,2);
        
        cnt = 0;
        for dx = (-radius):radius
            for dy = (-radius:radius)
                xx = x + dx;
                yy = y + dy;
                if (xx >= 1 && xx <= gmm.nX) && (yy >= 1 && yy <= gmm.nY)
                    % X_ are stored row-major.
                    % TODO: Perhaps revise?
                    iy = gmm.nY * (xx - 1) + yy;
                    if gmm.s_z(iy) == 1 % background
                        cnt = cnt + 1;
                    end
                end
            end
        end
        if cnt == 0
            cnt = cnt + 1; % smooth
        end
        like(n) = cnt / mass;    
    end
end
