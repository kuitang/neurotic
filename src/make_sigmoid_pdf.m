function [ ss ] = make_sigmoid_pdf( center, precision )
    pre_ss = @(x) 1-(1./(1+exp(-precision * (x - center))));
    Z = quad(pre_ss, 0, 1);
    ss = @(x) 1/Z * pre_ss(x);
end

