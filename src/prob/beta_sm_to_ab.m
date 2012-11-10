function [ a b ] = beta_sm_to_ab( s, m)
    a = s .* m;
    b = s .* (1 - m);
end

