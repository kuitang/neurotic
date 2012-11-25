function [ s m ] = beta_ab_to_sm( a, b )
% Convert standard beta parameterization to scale-mean
    s = a + b;
    m = a ./ s;
end

