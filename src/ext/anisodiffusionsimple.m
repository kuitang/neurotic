function anisodiffusionsimple(a)
a = im2double(a);
g1 =a;
g2 =a;
stepsize1=.2;
stepsize2=.1;
nosteps=200;
verbose=1;
if verbose
    figure(verbose);
    subplot(2,3,1); imshow(g1); title('Original Image'); drawnow;
end
h=waitbar(0,'Now is doing Diffusion');
for i=1:nosteps
    laplaceg1=snldStep1(g1);
    g1 = g1 +stepsize1 * laplaceg1;
    laplaceg2=snldStep2(g2);
    g2 = g2 +stepsize2 * laplaceg2;
    waitbar(i/nosteps,h);
    if verbose
        figure(verbose);
        subplot(2,3,2); imshow(uint8(255*(laplaceg1-min(laplaceg1(:)))/(max(laplaceg1(:))-min(laplaceg1(:)))));
        %subplot(2,2,2); imshow(uint8(255*(laplaceg-min(laplaceg(:)))/(max(laplaceg(:))-min(laplaceg(:)))));%????
        title('Laplace term')
        subplot(2,3,3); imshow(g1);
        title('Anisotropic Diffusion'); drawnow;
        subplot(2,3,5); imshow(uint8(255*(laplaceg2-min(laplaceg2(:)))/(max(laplaceg2(:))-min(laplaceg2(:)))));
        %subplot(2,2,2); imshow(uint8(255*(laplaceg-min(laplaceg(:)))/(max(laplaceg(:))-min(laplaceg(:)))));%????
        title('Laplace term')
        subplot(2,3,6); imshow(g2);
        title('Isotropic Diffusion'); drawnow;
    end
        waitbar(i/nosteps,h);
end
close(h);

function r = snldStep1( L )
% Discrete numerical scheme of dL/dt for scalar diffusion
N = size(L, 1);
M = size(L, 2);
lamda=.008;
% Set delta_x,delta_y=2,so we can set delta_x/2=delta_y/2=1
Lpc = translateImage( L, 1, 0 );
Lmc = translateImage( L, -1, 0 );
Lcp = translateImage( L, 0, 1 );
Lcm = translateImage( L, 0, -1 );
Lppc = translateImage( L, 2, 0 );
Lmmc = translateImage( L, -2, 0 );
Lcpp = translateImage( L, 0, 2 );
Lcmm = translateImage( L, 0, -2 );
Lap_i = (Lpc+Lcp)-2*L; % Laplace_I
d =abs( Lap_i/lamda );
e=-d.*d;
C = exp(e);
Cpc = translateImage( C, 1, 0 );
Cmc = translateImage( C, -1, 0 );
Ccp = translateImage( C, 0, 1 );
Ccm = translateImage( C, 0, -1);
r = ( 1/4*(Cpc.*( Lppc - L )-Cmc.*( L - Lmmc ))+1/4*(Ccp.*( Lcpp - L )-Ccm.*( L - Lcmm )) );

function k = snldStep2( L )
% Discrete numerical scheme of dL/dt for scalar diffusion
N = size(L, 1);
M = size(L, 2);
Lpc = translateImage( L, 1, 0 );
Lmc = translateImage( L, -1, 0 );
Lcp = translateImage( L, 0, 1 );
Lcm = translateImage( L, 0, -1 );
k = ( (Lpc-2*L+Lmc)+ (Lcp-2*L+Lcm) ); ...