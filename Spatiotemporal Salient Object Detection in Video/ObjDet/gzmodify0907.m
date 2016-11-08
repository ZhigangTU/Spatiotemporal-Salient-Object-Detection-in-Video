% Figure 5 in the paper
% robust batch image alignment example

% clear
clc ;
clear all;
close all ;

% addpath
addpath computeopticflow ;
addpath data ;
addpath results ;
%%
methodchoose =1, %original RPCA code
%methodchoose =5, 

%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;
pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates
userName = 'smtreeman100320';

% output path
if methodchoose==1
    destRoot = fullfile(currentPath,'orirpcaresults') ;
end   

if methodchoose==6
    destRoot = fullfile(currentPath,'gzblock2and1norm') ;
end

destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;
%% read every images
testImage = imread(fileNames{1});

if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);
%for color image
%w=w/3;

data = zeros([h*w, numImages], 'double');
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  if isrgb(currentImage)
    currentImage = rgb2gray(currentImage);
  end
  im = reshape(currentImage,w*h,1);
  data(:, fileIndex) = im(:,1);
end
data=data/255.0;

if methodchoose==1
    [A_hat E_hat iter] = inexact_alm_rpca(data);%original ALM method
    %[A_hat E_hat iter] = inexact_alm_rpca(data, 1 / sqrt(160*120));%original ALM method
    
    %to use the relative chage to obtain masks
%     RelativeChange = gzComputeRelativeChange(data, A_hat);
%     [MaskBlock1Pixle2 E_blockmask]  = gzEmatrixFindMasksChangingblocksize(RelativeChange, 0.2, 4);
    
    %to use E matrix to obtain masks
    Eabs = abs(E_hat);
    [MaskBlock1Pixle2 E_blockmask]  = gzEmatrixFindMasksChangingblocksize(Eabs, 0.2, 4);
end   

%[Mask_matrix E_blockmask]=gzAnalyzeEmatrix(E_hat);
%blocksize = 4 case


%%save the data as images (low rank and sparse matrix)
for fileIndex = 1:numImages
    immask = Mask_matrix(:,fileIndex);
    immask = reshape(immask,h,w);
    immask = immask * 255;
    %not display now
    %figure;
    %imshow(uint8(im1));
    outputMask  = sprintf('mask%d.jpg',fileIndex);
    outputMasks = fullfile(destDir, outputMask); 
    imwrite(immask, outputMasks);
    
    
    im1 = E_hat(:,fileIndex);
    for j = 1:h*w,
%         if (abs(im1(j,1))) < 0.02,
% 	        im1(j,1) = 0;
%         else
%             im1(j,1)=255;
%         end
        
        if im1(j,1) < 0
	      im1(j,1) = 255;
        end
    end
    im1 = reshape(im1,h,w);
    %not display now
    %figure;
    %imshow(uint8(im1));
    
    outputFileNameLowrank  = sprintf('Lowrank%d.jpg',fileIndex);
    outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
    imwrite(im1, outputFileNamesLowrank);
    
    
    im2 = A_hat(:,fileIndex);
    im2 = reshape(im2,h,w);
    im2 = im2 * 255;
    %not display now
    %figure;
    %imshow(uint8(im2));
    
    outputFileNameSparse  = sprintf('Sparse%d.jpg',fileIndex);
    outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
    %imwrite(im2, outputFileNamesSparse);
    imwrite(uint8(im2), outputFileNamesSparse);
end
