function [ h ] = visualize_gamma( shapes, scales )

    N = length(shapes);
    M = length(scales);
    
    series = cell(1, 2*N*M);
    leg    = cell(1, N*M);
    
    i = 1;
    for n = 1:N
        for m =1:M
            shape = shapes(n);
            rate  = 1 / scales(m);
            stop = gaminv(0.999, shape, rate);
            X = linspace(0, stop);
            Y = gampdf(X, shape, rate);
            series{2*i-1} = X;
            series{2*i}     = Y;
            leg{i} = [ 'shape = ' num2str(shape) '; scale = ' num2str(scales(m)) ];
            i = i + 1;
        end
    end
    
    figure
    hold on
    title('Behavior of Gamma');
    plot(series{:});
    legend(leg{:});        

end

