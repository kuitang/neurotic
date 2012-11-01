function [ samples ] = run_mm(iter_f, X, img, mm, niters, burn_in, sample_freq, recovery_file )    
    
    [N, D] = size(X);
    samples = [];
    mm.n = accumarray(mm.s_z, ones(1,N));    
    
    if nargin >=8 && ~isempty(recovery_file)
        load(recovery_file);
        first = iter + 1;        
    else
        first = 1;
    end
        
    for i = first:niters
        tic
        z_old = mm.s_z;
                
        mm = iter_f(mm, X);
        change = mean(abs(z_old ~= mm.s_z));
        
        itertime = toc;
        rate     = N / itertime;        
        
        fprintf(1, 'Iter %d: time = %0.2f; speed = %.2f; change = %0.3f; loglike = %06.2f\n', ...
                i, itertime, rate, change, mm.loglike);     
        
        if mod(i, sample_freq) == 0
            plot_point_overlay(img, mm, X);
            drawnow;
            
            if i > burn_in
                samples = [samples mm];
                iter = i;

                %save(['iter' num2str(i) '.mat']);
            end
        end        

    end
    
end


