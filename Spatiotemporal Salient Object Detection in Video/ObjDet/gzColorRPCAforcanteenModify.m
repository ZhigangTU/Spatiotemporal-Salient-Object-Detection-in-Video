% Gao Zhi rpca for color sequence, to obatin background
% robust batch image alignment example

% clear
clc ; clear all; close all ;

% addpath
addpath data ;
addpath results ;
%%
methodchoose =1001, %original RPCA code
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;

%userName = 'artscanteen0418ready';
%userName = 'manyBirds';
%userName = '5ready';
%userName = '2ready';
%userName = '0426morningready';
%userName = '0426morningreadyfortest';
%userName = '1testready';

%userName = '0503afterlunchpeakready';
%userName = '0503beforelunchpeakready';
%userName = '0503right2closeready';

%userName = '0502rainyafterlunchpeakready';
%userName = '0502rainybeforelunchpeakready';
%userName = '0502rainyafterlunchpeakagainready';

%userName = '04240502rainyafterlunchpeakagainready';
%userName = '0502rainybeforelunchpeakready';
userName = 'temp4';








% output path
if methodchoose==1001
    destRoot = fullfile(currentPath,'RPCAresultsoriginalpart') ;
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

dataR = zeros([h*w, numImages], 'double');
dataG = zeros([h*w, numImages], 'double');
dataB = zeros([h*w, numImages], 'double');
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  
  currentImageR = currentImage(:,:,1);
  currentImageG = currentImage(:,:,2);
  currentImageB = currentImage(:,:,3);
  
%   if isrgb(currentImage)
%     currentImage = rgb2gray(currentImage);
%   end
  
  imR = reshape(currentImageR,w*h,1);
  dataR(:, fileIndex) = imR(:,1);
  
  imG = reshape(currentImageG,w*h,1);
  dataG(:, fileIndex) = imG(:,1);
  
  imB = reshape(currentImageB,w*h,1);
  dataB(:, fileIndex) = imB(:,1);
end

dataR=dataR/255.0;
dataG=dataG/255.0;
dataB=dataB/255.0;

if methodchoose==1001
    [A_hatR E_hatR iterR] = inexact_alm_rpca(dataR);%original ALM method
    %clear E_hatR;
    clear iterR;
    %clear A_hatR
    %clear dataR;
    
    [A_hatG E_hatG iterG] = inexact_alm_rpca(dataG);%original ALM method
    %clear E_hatG;
    clear iterG;
    %clear A_hatG
    %clear dataG;
    
    [A_hatB E_hatB iterB] = inexact_alm_rpca(dataB);%original ALM method
    %clear E_hatB;
    clear iterB;
    %clear A_hatB
    %clear dataB;
    
end   

dRatio = 1.0;%%0.5;%%
E_hatR = abs(E_hatR);
meanR  = mean2(E_hatR);
stdR   = std2(E_hatR);
ThreshR=meanR+dRatio*stdR;
E_hatR=(E_hatR<ThreshR);

E_hatG = abs(E_hatG);
meanG  = mean2(E_hatG);
stdG   = std2(E_hatG);
ThreshG=meanG+dRatio*stdG;
E_hatG=(E_hatG<ThreshG);

E_hatB = abs(E_hatB);
meanB  = mean2(E_hatB);
stdB   = std2(E_hatB);
ThreshB=meanB+dRatio*stdB;
E_hatB=(E_hatB<ThreshB);

MaskSafe=E_hatR.*E_hatG;
MaskSafe=MaskSafe.*E_hatB;

MaskSafeAnti= 1-MaskSafe;


dataR =dataR.*MaskSafe+MaskSafeAnti;
dataG =dataG.*MaskSafe+MaskSafeAnti;
dataB =dataB.*MaskSafe+MaskSafeAnti;

    
% AANames = fullfile(destDir, 'Amatrix.mat'); 
% save(AANames,'A_hat') ;
% EENames = fullfile(destDir, 'Ematrix.mat'); 
% save(EENames,'E_hat') ;

%%[Mask_matrix E_blockmask]=gzAnalyzeEmatrix(E_hat);%%gaozhi

%%save the data as images (low rank and sparse matrix)
for fileIndex = 1:numImages
%     immask = Mask_matrix(:,fileIndex);
%     immask = reshape(immask,h,w);
%     immask = immask * 255;
%     %not display now
%     %figure;
%     %imshow(uint8(im1));
%     outputMask  = sprintf('mask%d.jpg',fileIndex);
%     outputMasks = fullfile(destDir, outputMask); 
%     imwrite(immask, outputMasks);
    


    imRbackground = A_hatR(:,fileIndex);
    imRbackground = reshape(imRbackground,h,w);
    imRbackground = imRbackground * 255;
    imR = dataR(:,fileIndex);
    imR = reshape(imR,h,w);
    imR = imR * 255;
    
    imGbackground = A_hatG(:,fileIndex);
    imGbackground = reshape(imGbackground,h,w);
    imGbackground = imGbackground * 255;
    imG = dataG(:,fileIndex);
    imG = reshape(imG,h,w);
    imG = imG * 255;
    
    imBbackground = A_hatB(:,fileIndex);
    imBbackground = reshape(imBbackground,h,w);
    imBbackground = imBbackground * 255;
    imB = dataB(:,fileIndex);
    imB = reshape(imB,h,w);
    imB = imB * 255;
    
    BkImage(:,:,1)=imR;
    BkImage(:,:,2)=imG;
    BkImage(:,:,3)=imB;
    LowrankImage(:,:,1)=imRbackground;
    LowrankImage(:,:,2)=imGbackground;
    LowrankImage(:,:,3)=imBbackground;
    
    outputFileNameSparse  = sprintf('Trainingarea%d.jpg',fileIndex);
    outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
    %imwrite(im2, outputFileNamesSparse);
    %imwrite(uint8(im2), outputFileNamesSparse);
    imwrite(uint8(BkImage), outputFileNamesSparse);
    
    
    outputFileNameLowrank  = sprintf('Lowrank%d.jpg',fileIndex);
    outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
    %imwrite(im2, outputFileNamesSparse);
    %imwrite(uint8(im2), outputFileNamesSparse);
    imwrite(uint8(LowrankImage), outputFileNamesLowrank);
end
