% clear
clc ;
clear all;
close all ;

% addpath
addpath data ;
addpath results ;
%% define images' path

currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;
userName = 'smtreeman100320' ;

% output path
destRoot1 = fullfile(currentPath,'split4subs') ;
destRoot = fullfile(destRoot1,userName) ;

destDir1 = fullfile(destRoot,'sub1') ;
if ~exist(destDir1,'dir')
    mkdir(destRoot,'sub1') ;
end

destDir2 = fullfile(destRoot,'sub2') ;
if ~exist(destDir2,'dir')
    mkdir(destRoot,'sub2') ;
end

destDir3 = fullfile(destRoot,'sub3') ;
if ~exist(destDir3,'dir')
    mkdir(destRoot,'sub3') ;
end

destDir4 = fullfile(destRoot,'sub4') ;
if ~exist(destDir4,'dir')
    mkdir(destRoot,'sub4') ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;
%% read every images
testImage = imread(fileNames{1});

if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);
   
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  if isrgb(currentImage)
    currentImage = rgb2gray(currentImage);
  end
  
  Submatrix11  = currentImage(1:h/2, 1:w/2);
  outputImgUL = sprintf('UL%05d.jpg',fileIndex);
  outputUL    = fullfile(destDir1, outputImgUL); 
  imwrite(Submatrix11, outputUL);
  
  Submatrix12  = currentImage(1:h/2, 1+w/2:w);
  outputImgUR = sprintf('UR%05d.jpg',fileIndex);
  outputUR    = fullfile(destDir2, outputImgUR); 
  imwrite(Submatrix12, outputUR);
  
  Submatrix21  = currentImage(1+h/2:h, 1:w/2);
  outputImgDL = sprintf('DL%05d.jpg',fileIndex);
  outputDL    = fullfile(destDir3, outputImgDL); 
  imwrite(Submatrix21, outputDL);
  
  Submatrix22  = currentImage(1+h/2:h, 1+w/2:w);
  outputImgDR = sprintf('DR%05d.jpg',fileIndex);
  outputDR    = fullfile(destDir4, outputImgDR); 
  imwrite(Submatrix22, outputDR);
end