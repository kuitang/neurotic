function [ a b ] = beta_sm_to_ab( s, m)
% Convert scale-mean parameterization of Beta distribution to standard
    a = s .* m;
    b = s .* (1 - m);
end

