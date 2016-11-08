% Gao Zhi read canteen (or other) color image sequence, 
% and extract useful part, and save it for rpca 

% clear
clc ; clear all; close all ;

% addpath
% addpath dataNeedPrepare ;
% addpath results ;
addpath thermal ;

%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'thermal') ;

%userName = 'traffic';
%userName = 'sidewalk';
%userName = 'original00900onwards';
%userName = 'corridor';
%userName = 'canoeinput';%park

%userName = 'original00900onwards';
% userName = 'original00900onwardsparts';
%userName = 'fountain01input';
% userName = 'parkthermal240by320';
userName = 'canoe240by320';




% output path
destRoot = fullfile(currentPath,'thermalresize') ;
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

%%
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  %testImage = currentImage(1:h-20,:,:);
  %ReSizeIMG = imresize(testImage, [240 320]);
  %ReSizeIMG = imresize(testImage, [60 80]);
  %currentImage = currentImage(80:h,1:320,:);
  %ReSizeIMG = imresize(currentImage, [60 80]);
  %ReSizeIMG = currentImage(150:209,8:87,:);
  
  %ReSizeIMG = imresize(currentImage, [240 320]);
  %ReSizeIMG = imresize(currentImage, [60 80]);
  
  ReSizeIMG = imresize(currentImage, [120 160]);
  
  partnames = sprintf('GZ%06d.jpg',fileIndex);
  outputFileName  = strcat(userName, partnames); 
  %outputFileName  = sprintf('GZtraffic%06d.jpg',fileIndex);
  outputFileNames = fullfile(destDir, outputFileName); 
  %imwrite(im2, outputFileNamesSparse);
  %imwrite(uint8(im2), outputFileNamesSparse);
  %imwrite(uint8(BkImage), outputFileNamesSparse);
  imwrite(ReSizeIMG, outputFileNames);
end

