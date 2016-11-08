%reorder the file, and rename it.

% clear
clc ; clear all; close all ;

addpath data ;
addpath datarenameoutput ; 

currentPath = cd;
% input path
imagePath = fullfile(currentPath,'data') ;

%userName  = 'cam2new';
userName  = 'cam1on20130125';
%userName  = 'cam1on13July';

destRoot = fullfile(currentPath,'datarenameoutput');
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  
  %currentImage = imresize(currentImage, [240 320]);
 
  outputFileName  = sprintf('ColorPrepare%06d.jpg',fileIndex);
  outputFileNames = fullfile(destDir, outputFileName); 
  
  imwrite(currentImage, outputFileNames);
end