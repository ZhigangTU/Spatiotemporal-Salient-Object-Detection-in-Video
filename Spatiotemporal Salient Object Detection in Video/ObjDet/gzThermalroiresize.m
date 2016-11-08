% Gao Zhi read canteen (or other) color image sequence, 
% and extract useful part, and save it for rpca 

% clear
clc ; clear all; close all ;

% addpath
% addpath dataNeedPrepare ;
% addpath results ;
addpath thermal ;
addpath thermalroi ;


%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'thermal') ;
imagePathroi = fullfile(currentPath,'thermalroi') ;

%userName = 'traffic';
%userName = 'sidewalk';
%userName = 'original00900onwards';
%userName = 'corridor';
%userName = 'canoeinput';%park

%userName = 'original00900onwards';
% userName = 'original00900onwardsparts';

% userName = 'streetLightinput';
%userName = 'streetlightgroundtruth175';
userName = 'streetlightgroundtruth175roi';





% output path
% destRoot = fullfile(currentPath,'thermalresize') ;
destRoot = fullfile(currentPath,'thermalroiresize') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

[fileNamesroi, numImagesroi] = gzget_training_images( imagePathroi, userName) ;

%% read msk image / just one image
maskImage = imread(fileNamesroi{1});
if isrgb(maskImage)
    maskImage = rgb2gray(maskImage);
end
maskBinary = im2bw(maskImage, 0.5);
maskBinary = im2bw(maskImage, 0.5)*1.0;
% imshow(maskBinary);
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
%for fileIndex = 175:3:numImages
  currentImage = imread(fileNames{fileIndex});
  
  if isrgb(currentImage)
    currentImage = rgb2gray(currentImage);
  end
  %currentImage=currentImage*1.0;
  currentImage=double (currentImage);
  currentImage=currentImage.*maskBinary;
  
  %imshow(currentImage/255);
  
  %testImage = currentImage(1:h-20,:,:);
  %ReSizeIMG = imresize(testImage, [240 320]);
  %ReSizeIMG = imresize(testImage, [60 80]);
  currentImage = currentImage(30:170,230:320);
  ReSizeIMG = imresize(currentImage, [60 80]);
  %ReSizeIMG = currentImage(150:209,8:87,:);
  
  %ReSizeIMG = imresize(currentImage, [240 320]);
  %ReSizeIMG = imresize(currentImage, [60 80]);
  
  %partnames = sprintf('GZ%06d.jpg',fileIndex);
  partnames = sprintf('GZ%06d.png',fileIndex);
  outputFileName  = strcat(userName, partnames); 
  %outputFileName  = sprintf('GZtraffic%06d.jpg',fileIndex);
  outputFileNames = fullfile(destDir, outputFileName); 
  %imwrite(im2, outputFileNamesSparse);
  %imwrite(uint8(im2), outputFileNamesSparse);
  %imwrite(uint8(BkImage), outputFileNamesSparse);
  imwrite(uint8(ReSizeIMG), outputFileNames);
  
end

