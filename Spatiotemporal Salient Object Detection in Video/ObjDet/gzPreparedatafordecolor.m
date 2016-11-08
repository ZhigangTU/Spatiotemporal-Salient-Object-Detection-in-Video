% gaozhi read img sequence and save .mat file for DECOLOR processing
% clear
clc ;
clear all;
close all ;

% addpath
addpath data ;
addpath results ;
%%
%methodchoose =1, %original RPCA code
%methodchoose   = 21,%block RPCA method
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'gz img sequence') ;
%userName = 'rainca60' ;
%userName = 'treeman50';
%userName = 'cardimlight60';
%userName = 'scienceblue';
%userName = 'basicpart';
%userName = 'lightswitchshanmo';
%userName = 'scienceBigRain';
%userName = 'flock100';
% userName = 'boats';
% userName = 'sidewalkallspecialpart2';
% userName = 'aligned00900aligned001to294range00085to00157';
%userName = 'aligned00900aligned295after00120to00210';
% userName = 'originalsidewalkallspecialpart2';
% userName = 'original00900onwardsaligned001to294from85to157';
% userName = 'original00900onwardsaligned295after120to210';
% userName = 'parkthermal240by320from440to540';
% userName = 'canoe240by320from835to970';
userName = 'canoe240by320from750to1000';





%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

%% read first image to get the size info
testImage = imread(fileNames{1});
if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

%% read image file, to establish the data matrix
hResize = 120;
wResize = 160;
ImData = zeros([hResize, wResize, numImages], 'double');
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  if isrgb(currentImage)
    currentImage = rgb2gray(currentImage);
  end
  %im = reshape(currentImage,w*h,1);
  % gaozhi add the resize step
  currentImage = imresize(currentImage, [hResize wResize]);
  ImData(:,:, fileIndex) = currentImage;
  
  %outputtmp  = sprintf('img%05d.jpg',fileIndex);
  %imwrite(currentImage, outputtmp);
  
end

ImData=ImData/255.0;


    

% output path

destDir = fullfile(imagePath,userName) ;
FileNames = strcat(userName,'.mat'); 
SaveNames = fullfile(destDir, FileNames); 
save(SaveNames,'ImData') ;