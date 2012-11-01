function [ y ] = log1pexp( x )
% Generalized log1p y \approx log(1 + exp(x)). Uses y = log1p(x) for small
% x and y = x for large x.

    small_idxs = x < 100;
    y = zeros(size(x));
    y(small_idxs)  = log1p(exp(x(small_idxs)));
    y(~small_idxs) = x(~small_idxs);    

end

