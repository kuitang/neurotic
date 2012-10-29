function [ p ] = fast_mvtpdf( mvtparams, x )
    % We scale by both mean and the correlation magnitude, which MATLAB
    % does not.
    dx = x - mvtparams.mean;    
    dx = dx * mvtparams.stdev'; % MATLAB's chol is L'L!
    
    p = mvtparams.Z * (1 + sum(dx.^2) / mvtparams.dof) ^ mvtparams.pow;
end

