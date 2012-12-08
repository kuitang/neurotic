%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getEdgeFeatures.m picks edge enhancing Gaussian Second Derivative 
%                   Filters
%
% Copyright (c) Ritwik Kumar, Harvard University 2010
%               www.seas.harvard.edu/~rkkumar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KT - This is R(x,y) (equation 4) in the paper. Convolution becomes
% multiplication in frequency domain.

function [sf] = getEdgeFeatures(im, sup, sc, or)

[imr imc] = size(im);
imi = fft2(double(im));
S=makeBarFilters(sup,or,sc);
nof = size(S,3);
sf = zeros(imr,imc,nof,'single');
for i=1:nof
    sf(:,:,i)= (circshift(ifft2(imi.*fft2(S(:,:,i),imr,imc)),[-fix(sup/2) -fix(sup/2)]));
end


