% Gao Zhi read canteen color image sequence, and extract useful part, and save it for rpca 

% clear
clc ; clear all; close all ;
% addpath
%addpath data ;
%addpath results ;
addpath wave ;
%%
methodchoose =101, %original RPCA code
%% define images' path
currentPath = cd;

% input path
% imagePath = fullfile(currentPath,'data') ;
imagePath = fullfile(currentPath,'wave') ;

%userName = '0502rainyafterlunchpeak';
% userName = 'WhimbrelsSungeiBulohWetlandReserves';
% userName = 'AnimalCamCrocodile';
%userName = 'Frontier_Cam02good130120';
% userName = 'CivetCat';
% userName = 'AnimalCamJumpingfish';
% userName = 'GreyHeronSungeiBulohWetlandReserve';
% userName = 'GiantMudskipper';
% userName = 'YellowBittern';
% userName = '1HourRealNatural09363to12485';
userName = 'DS_0068';

nSampleJump=7;


% output path
if methodchoose==101
%     destRoot = fullfile(currentPath,'colorprepareforrpca') ;
destRoot = fullfile(currentPath,'wavegoodpart') ;
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

%dataR = zeros([h*w, numImages], 'double');
%dataG = zeros([h*w, numImages], 'double');
%dataB = zeros([h*w, numImages], 'double');
nCount =1;

for fileIndex = 1:nSampleJump:numImages
    currentImage = imread(fileNames{fileIndex});
    testImage = currentImage(5:h-5,10:w-10,:);
    % testImage = currentImage(:,:,:);
    ReSizeIMG = imresize(testImage, [240 320]);
    %   ReSizeIMG = imresize(testImage, [60 80]);

    %outputFileName  = sprintf('ColorPrepare%06d.jpg',fileIndex);
    %   outputFileName  = sprintf('Training%06d.jpg',fileIndex);
    %   outputFileName  = sprintf('add%06d.jpg',fileIndex);
%     outputFileName  = sprintf('in%06d.jpg',fileIndex);
    outputFileName  = sprintf('in%06d.jpg',nCount);
    outputFileNames = fullfile(destDir, outputFileName);
    %imwrite(im2, outputFileNamesSparse);
    %imwrite(uint8(im2), outputFileNamesSparse);
    %imwrite(uint8(BkImage), outputFileNamesSparse);
    imwrite(ReSizeIMG, outputFileNames);
    
    nCount=nCount+1;
  
end
