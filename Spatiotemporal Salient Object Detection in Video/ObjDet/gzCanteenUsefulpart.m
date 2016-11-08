% Gao Zhi read canteen (or other) color image sequence, 
% and extract useful part, and save it for rpca 

% clear
clc ; clear all; close all ;

% addpath
addpath dataNeedPrepare ;
%addpath results ;

%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'dataNeedPrepare') ;

%userName = '0506sunday';
%userName = '0502rainydaypartofall';
%userName = 'traffic';
%userName = 'sidewalk';
%userName = 'original00900onwards';
userName = 'sidewalkall';






% output path
destRoot = fullfile(currentPath,'GetColorSequenceforRPCA') ;
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
  testImage = currentImage(1:h-20,:,:);
  %ReSizeIMG = imresize(testImage, [240 320]);
  ReSizeIMG = imresize(testImage, [60 80]);
  
  %ReSizeIMG = imresize(currentImage, [240 320]);
  %ReSizeIMG = imresize(currentImage, [60 80]);
  
  partnames = sprintf('GZ%06d.jpg',fileIndex);
  outputFileName  = strcat(userName, partnames); 
  %outputFileName  = sprintf('GZtraffic%06d.jpg',fileIndex);
  outputFileNames = fullfile(destDir, outputFileName); 
  %imwrite(im2, outputFileNamesSparse);
  %imwrite(uint8(im2), outputFileNamesSparse);
  %imwrite(uint8(BkImage), outputFileNamesSparse);
  imwrite(ReSizeIMG, outputFileNames);
  
end





%% read every images
% testImage = imread('SceneSampling000001.jpg');
% 
% if isrgb(testImage)
%     grayImage = rgb2gray(testImage);
% end
% 
% [h, w] = size(testImage);
% 
% testImage = testImage(20:h,:,:);
% 
% ReSizeIMG = imresize(testImage, [240 320]);
% 
% imwrite(ReSizeIMG, 'Part000001.jpg');