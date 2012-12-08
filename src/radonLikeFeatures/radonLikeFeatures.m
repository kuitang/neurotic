function [ intx ] = radonLikeFeatures( imo, use_raw_imo, res, extract_func )
%radonLikeFeatures Generate radon-like features from scan segments
%
% intx = radonLikeFeatures(imo, res, extract_fun) given [X Y] image imo,
% returns [X Y res] feature responses in intx (one for each scan direction)
% 
% extract_func is a handle used in 
%   int(i,prev:j) = extract_func(im(i,prev:j));
%
% Adapted from radonLikeFeatureDemo.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radonLikeFeaturesDemo.m Demonstrates Radon-Like Features for the 
%             special case of edge enhancement.
%
% Copyright (c) Ritwik Kumar, Harvard University 2010
%               www.seas.harvard.edu/~rkkumar
%
% Please cite the following paper if you fine this code useful:
%  Ritwik Kumar, Amelio V. Reina & Hanspeter Pfister, “Radon-Like 
%  Features and their Application to Connectomics”, IEEE Computer 
%  Society Workshop on Mathematical Methods in Biomedical Image 
%  Analysis (MMBIA), 2010
%
% Various steps involved are explained in the comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear;clc;
% filename = 'sample.png'; % Input image
%res = 36; % # of directions of scanning
    ang = 360/res; % angular resolution of scanning
    %imo = double(imread(filename)); % read input

    %scaling the image to 0 - 255
    imo = imo - min(imo(:));
    imo = imo ./ max(imo(:)) * 255;
    [imro imco] = size(imo);

    %edge enhancing features
    %sf = getEdgeFeatures(imo,15,3,12);
    %imo = (max(sf,[],3));
    
    % In general, use the enhanced edge features to get knots only    
    sf = getEdgeFeatures(imo,15,3,12);
    knot_imo = (max(sf,[],3));
    
    if ~use_raw_imo
        imo = knot_imo;
    end

    %set the boundary to be 0 -- NO; SET IT TO BE WHITE!
    imo(1,:) = 0;imo(imro,:)=0;imo(:,1) = 0;imo(:,imco)=0;
    knot_imo(1,:) = 0;knot_imo(imro,:)=0;knot_imo(:,1) = 0;knot_imo(:,imco)=0;

    %set the knots
    edo = edge(knot_imo,'canny',0.37);    
    
    dumo = ones(imro,imco);% helper matrix to keep track of ROI

    %some variable initialization
    dum=1;ed=1;int=1;
    intx = zeros(imro,imco,res,1);

    for xx=1:res % repeat for each scan direction
        fprintf('Scan %d\n',xx);

        %rotate the image by the scanning angle
        %this allows us to scan always along the horizontal axis
        % NOT THE MOST EFFICIENT WAY OF DOING THING

        % Notation: the indices i,prev:j in the inenr loop below denotes a scan
        % segment (from knot to knot). (At least I think so because it looks
        % that way. Study the code in more detail.)

        im = imrotate(imo,ang*(xx-1)); 
        dum = imrotate(dumo,ang*(xx-1)); %rotate ROI tracker
        ed = imrotate(edo,ang*(xx-1)); %rotate knot function. 
        %This may change the know function slightly due to 
        %interpolation artifacts.

        [imr imc] = size(im);

        int = im .* 0;
        flag = 0;
        for i=1:imr %for each row in the image
            j=1;
            while((dum(i,j) == 0) && (j < imc)) %find ROI beginning
                j=j+1;
                flag = 1;
            end
            prev = j;
            while((dum(i,j) ~= 0) && (j < imc)) %till we get out of ROI
                if(ed(i,j) == 1) %We hit a know
                    if(j-prev < -1)
                        int(i,prev:j) = 255;
                    else
                        %extraction function which assigns mean value
                        %change this to change the extraction function
                        
                        int(i,prev:j) = extract_func(im(i,prev:j));
                        
                        %int(i,prev:j) = mean(im(i,prev:j)); 
                        %int(i,prev:j) = min(im(i,prev:j));
                    end
                    prev = j;
                    flag = 0;
                end
                j=j+1;
            end
        end

        int = imrotate(int,-ang*(xx-1)); %rotate the feature image back
        [imr imc] = size(int);
        osetr = floor((imr-imro)/2);
        osetc = floor((imc-imco)/2);
        int = int(osetr+1:osetr+imro,osetc+1:osetc+imco); %extract ROI
        intx(:,:,xx) = int; %store ROI
    end
    %show the mean of the computed Radon-Like Features
    
    %imshow(mean(intx,3),[]);


end

