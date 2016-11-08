% Gao Zhi read canteen color image sequence, and extract useful part, and save it for rpca 

% clear
clc ; clear all; close all ;

% addpath
addpath data ;
%addpath results ;
%%
methodchoose =201, %original RPCA code
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;

%userName = 'artscanteen0418';
%userName = 'manyBirds';
%userName = '2';
%userName = 'toreorderbackground';
%userName = '0502merge0424'


userName = '0502dataallready';


% output path
if methodchoose==201
    destRoot = fullfile(currentPath,'Reorderfile') ;
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

for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  %testImage = currentImage(20:h,:,:);
  %ReSizeIMG = imresize(testImage, [240 320]);
  
  outputFileName  = sprintf('backgroundtraining%06d.jpg',fileIndex+312);
  %outputFileName  = sprintf('backgroundtraining%06d.jpg',fileIndex);
  outputFileNames = fullfile(destDir, outputFileName); 
  %imwrite(im2, outputFileNamesSparse);
  %imwrite(uint8(im2), outputFileNamesSparse);
  %imwrite(uint8(BkImage), outputFileNamesSparse);
  imwrite(currentImage, outputFileNames);
  
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