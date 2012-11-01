function [ fig, match_idx ] = plot_sift( I1, I2, d1, d2, f1, f2, matches, scores )

[~, all_match_idx] = sort(scores);
    
    if length(all_match_idx > 100)
        match_idx = all_match_idx(1:100);
    else
        match_idx = all_match_idx;
    end
    N = length(match_idx);
    c = jet(N);

    fig = figure('Units', 'Normalized', 'OuterPosition',[0 0.3 1 0.7])
    title(['bock11_begin, slices 1 and 2, ', num2str(length(scores)) ' matched SIFT points']);

    subplot(121);
    imagesc(I1);
    colormap(gray);
    for i = 1 : N
        h = vl_plotsiftdescriptor(d1(:,matches(1,match_idx(i))), f1(:,matches(1,match_idx(i))));
        set(h, 'color', c(i,:));    
    end

    subplot(122);
    imagesc(I2);
    colormap(gray);
    for i = 1 : N
        h = vl_plotsiftdescriptor(d2(:,matches(2,match_idx(i))), f2(:,matches(2,match_idx(i))));    
        set(h, 'color', c(i,:));    
    end

end

