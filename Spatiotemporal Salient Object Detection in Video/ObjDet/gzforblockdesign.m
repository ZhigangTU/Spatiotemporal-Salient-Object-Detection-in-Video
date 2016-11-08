% Figure 5 in the paper
% robust batch image alignment example

% clear
clc ;
clear all;
close all ;

% addpath
addpath data ;
addpath results ;
%%
%methodchoose =1, %original RPCA code
%methodchoose =5, 
%             =2, try to enforce Block constraint
%methodchoose  =3;% basic mix 2 and 1 norm regularizer
%methodchoose  =4,
%methodchoose   =6,% block 2 1 norm 
%methodchoose   =7,
methodchoose   =8,
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;
pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates
userName = 'smtreeman100320' ;

% output path
if methodchoose==8
    destRoot = fullfile(currentPath,'gztodesignblock') ;
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

if methodchoose==8
    [A_hat E_hat iter] = gzblock2and1norm_inexact_alm_rpca(data);%
    [Mask_matrix E_blockmask]=gzAnalyzeEmatrix(E_hat);
    %[A_hat E_hat iter] = gzAdaptiveThreshblock2and1norm_inexact_alm_rpca(data, E_blockmask);%
end 

AANames = fullfile(destDir, 'Amatrix.mat'); 
save(AANames,'A_hat') ;
EENames = fullfile(destDir, 'Ematrix.mat'); 
save(EENames,'E_hat') ;

%%save the data as images (low rank and sparse matrix)
for fileIndex = 1:numImages
    immask = Mask_matrix(:,fileIndex);
    immask = reshape(immask,h,w);
    immask = immask * 255;
    %not display now %figure; %imshow(uint8(im1));
    outputMask  = sprintf('mask%d.jpg',fileIndex);
    outputMasks = fullfile(destDir, outputMask); 
    imwrite(immask, outputMasks);
end
