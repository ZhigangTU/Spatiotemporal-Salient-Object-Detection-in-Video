% Gao Zhi read canteen color image sequence, and extract useful part, and save it for rpca 

% clear
clc ; clear all; close all ;
% addpath
addpath data ;
%addpath results ;
%%
methodchoose =101, %original RPCA code
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;

%userName = '0502rainyafterlunchpeak';
userName = 'trafficpart120to210';


% output path
if methodchoose==101
    destRoot = fullfile(currentPath,'forRASL') ;
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
  %testImage = currentImage(1:224,:,:);
  %testImage = currentImage(:,:,:);
  testImage = rgb2gray(currentImage);
  %ReSizeIMG = imresize(testImage, [120 160]);
  ReSizeIMG = imresize(testImage, [60 80]);
  
  %outputFileName  = sprintf('ColorPrepare%06d.jpg',fileIndex);
%   outputFileName  = sprintf('Training%06d.jpg',fileIndex);
  outputFileName  = sprintf('In%06d.jpg',fileIndex);
  outputFileNames = fullfile(destDir, outputFileName); 
  %imwrite(im2, outputFileNamesSparse);
  %imwrite(uint8(im2), outputFileNamesSparse);
  %imwrite(uint8(BkImage), outputFileNamesSparse);
  imwrite(ReSizeIMG, outputFileNames);
  
end

