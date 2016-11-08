% gaozhi wang yuxiang modified, for RASL robust batch image alignment example
 
clc; 
clear all; 
close all ;

% addpath
addpath RASL_toolbox ;
addpath data ;
addpath results ;

%% define images' path
currentPath = cd;
% input path
imagePath = fullfile(currentPath,'data');
pointPath = fullfile(currentPath,'data');  % path to files containing initial feature coordinates

% userName = 'carheavyrain';
% userName = 'gzAl_Gore';
% userName = 'Al_Gore';
% userName = 'traffic80by60num300';
% userName = 'traffic80by60num300eccvback';
% userName = 'trafficpart120to210with60by80';
% userName = 'DS_0049original60by80';
% userName = 'DS_0004original60by80';
userName = 'DS_0004original60by80';

% output path
destRoot = fullfile(currentPath,'results') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

% datagz = zeros([60, 80], 'double'); imshow(datagz);

%% define parameters

% display flag
raslpara.DISPLAY = 1 ;

% save flag
raslpara.saveStart = 1 ;
raslpara.saveEnd = 1 ;
raslpara.saveIntermedia = 0 ;

% for face images
%raslpara.canonicalImageSize = [ 80 60  ];
raslpara.canonicalImageSize = zeros(1,2); % first for vertical second for horizontal
cropped_origin = zeros(1,2); % first for horizontal, second for vertical

cropped_origin = [5,5]; 
canonical_size = [50,70];
raslpara.canonicalImageSize =[50,70];
% raslpara.canonicalImageSize = canonical_imgsz/2;

% raslpara.canonicalCoords = [ 5  55 ; 32 32  ];
raslpara.canonicalCoords = [ 5  55 ; 32 32  ];
                            
% parametric tranformation model
% default 'AFFINE'; 
% raslpara.transformType = 'AFFINE'; 
% raslpara.transformType = 'HOMOGRAPHY';
% raslpara.transformType = 'TRANSLATION';
raslpara.transformType = 'EUCLIDEAN';
% one of 'TRANSLATION', 'EUCLIDEAN', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'

raslpara.numScales = 1 ; % if numScales > 1, we use multiscales

% main loop
raslpara.stoppingDelta = .01; % stopping condition of main loop
raslpara.maxIter = 25; % maximum iteration number of main loops

% inner loop
raslpara.inner_tol = 1e-6 ;
raslpara.inner_maxIter = 1000 ;
raslpara.continuationFlag = 1 ;
raslpara.mu = 1e-3 ;
raslpara.lambdac = 1.1 ; % lambda = lambdac/sqrt(m)

%% Get training images

% get initial transformation
% transformationInit = 'SIMILARITY';
% transformationInit = 'TRANSLATION';
transformationInit = 'IDENTITY';

[fileNames, transformations, numImages] = get_training_images( imagePath, pointPath, userName, raslpara.canonicalCoords, transformationInit) ;

transformations = {};
for i=1:numImages
    transformations{i} = eye(3)+ [zeros(3,2) [cropped_origin';0]];
end

%% RASL main loop: do robust batch image alignment
[D, Do, A, E, xi, numIterOuter, numIterInner ] = rasl_main(fileNames, transformations, numImages, raslpara, destDir);

%% plot the results
% layout.xI = 10 ;
% layout.yI = 14 ;
% layout.gap = 0 ;
% layout.gap2 = 0 ;
% rasl_plot(destDir, numImages, raslpara.canonicalImageSize, layout)

%% save the data as images (low rank and sparse matrix)
% h=80;
% w=60;
h=50;
w=70;
Dmax = max(max(D));
Domax = max(max(Do));
Amax =  max(max(A));

for fileIndex = 1:numImages
%     immask = Mask_matrix(:,fileIndex);
%     immask = reshape(immask,h,w);
%     immask = immask * 255;
%     % figure; imshow(uint8(im1));
%     outputMask  = sprintf('mask%d.jpg',fileIndex);
%     outputMasks = fullfile(destDir, outputMask); 
%     imwrite(immask, outputMasks);
    
    im1 = D(:,fileIndex)/Dmax;
    im1 = reshape(im1,h,w);
    im1 = im1 * 255;
    % not display now
    % figure; imshow(uint8(im1));
    
    outputFileNameLowrank  = sprintf('Input%06d.jpg',fileIndex);
    outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
    % imwrite(im1, outputFileNamesLowrank);
    imwrite(uint8(im1), outputFileNamesLowrank);
       
    im2 = Do(:,fileIndex)/Domax;
    im2 = reshape(im2,h,w);
    im2 = im2 * 255;
    outputFileNameSparse  = sprintf('Outputaligned%06d.jpg',fileIndex);
    outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
    %imwrite(im2, outputFileNamesSparse);
    imwrite(uint8(im2), outputFileNamesSparse);
    
    im3 = A(:,fileIndex)/Amax;
    im3 = reshape(im3,h,w);
    im3 = im3 * 255;
    outputFileNametmp  = sprintf('lowrankpart%06d.jpg',fileIndex);
    outputFileNamestmp= fullfile(destDir, outputFileNametmp); 
    %imwrite(im2, outputFileNamesSparse);
    imwrite(uint8(im3), outputFileNamestmp);
end
