function [ h ] = plot_thresh( I, intervals )

    edges = linspace(0, 1, intervals + 2);
    for i = 2:(length(edges) - 1)
        figure;
        imshow(I > edges(i));
        title(['threshold ' num2str(edges(i))]);
    end

end

