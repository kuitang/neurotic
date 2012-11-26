function test_online_cov(N, D, tol)
    % Notebook pp 41 (but we subtract a t)
    X = rand(N, D);

    mu = 0;
    Sigma = zeros(D);
    T = 0;

    % Addition: Update covariance (depends on current mean), then mean.
    %
    % The iteration begins at t-1 and ends at t in the equations.
    for n = 1:N    
        T = T + 1;
        dx = X(n,:) - mu;
        Delta = (T - 1)/T * (dx'*dx);
        Sigma = ((T - 1)*Sigma + Delta) / T;

        mu = mu + 1/T * dx;
    end

    assert(norm(mu - mean(X)) < tol);
    assert(norm(Sigma - cov(X, 1)) < tol);
    assert(T == N);

    % Subtraction: Update mean, then covariance (depends on "previous" mean,
    % hence "future" mean here.
    %
    % The iteration begins at t+1 and ends at t in the equations.
    for n = 2:N    
        mu = 1/(T - 1) * (T*mu - X(n,:));

        dx = X(n,:) - mu;
        Delta = (T - 1)/T * (dx'*dx);
        Sigma = (T*Sigma - Delta) / (T - 1);    

        T = T - 1;
    end    
    
    assert(T == 1);
    assert(norm(X(1,:) - mu) < tol);
    assert(norm(Sigma) < tol);
end
