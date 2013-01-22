function [ J ] = showslic( segments, I, mapzeroto )
    
    J = zeros(size(I), 'single');
    segments = segments + 1;
    N = max(segments(:));
    
    S = zeros(N, 1);
    n = zeros(N, 1);
    
    [X Y] = size(I);
    for j = 1:Y
        for i = 1:X
            seg = segments(i, j);
            S(seg) = S(seg) + I(i, j);
            n(seg) = n(seg) + 1;
        end
    end
    
    S = S ./ n;
    
    % For watershed instead of slic
    if nargin > 2
        S(1) = mapzeroto;
    end
    
    for j = 1:Y
        for i = 1:X                            
            J(i, j) = S(segments(i, j));
        end
    end    

end

