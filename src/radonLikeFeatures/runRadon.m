%% Configuration
res = 72;
img = f1;

%% Cell boundary enhancement

% Do the maximum Gaussian response
%sf = getEdgeFeatures(img,15,3,12);
%imo = (max(sf,[],3));

% argument 2 is use_raw_imo
boundary_intx = radonLikeFeatures(img, false, res, @(scanseg) mean(scanseg));
figure(1);
m_boundary_intx = mean(boundary_intx,3);
imshow(m_boundary_intx, []);

%% Cell background segmentation

bg_intx = radonLikeFeatures(img, true, res, @(scanseg) min(scanseg));
figure(2);
imshow(mean(bg_intx,3),[]);

% Plotting rescale
figure(3);
rescale_m = mean(bg_intx,3);
rescale_m = (rescale_m - min(rescale_m(:))) ./ (max(rescale_m(:)) - min(rescale_m(:)));
imshow(rescale_m);

%% Mitochondrial segmentation

mit_intx = radonLikeFeatures(img, true, res, @(scanseg) mean(scanseg));
figure(3);
imshow(mean(mit_intx,3),[]);
