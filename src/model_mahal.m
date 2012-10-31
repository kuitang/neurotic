function [ d ] = mahal( X, mu, S )
% Mahalobonis distance from model mean mu and covariance S.
    [N D] = size(X);
    dX = bsxfun(@minus, X, mu);
    [stdev, t] = cholcov(inv(S));
    assert(t == 0, 'S is not covariance!');
    dX = dX * stdev';
    d = sum(dX .* dX, 2);
end

