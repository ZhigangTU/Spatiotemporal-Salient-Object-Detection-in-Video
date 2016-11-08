% robust batch image alignment example

% clear
clc ; clear all; close all ;

% addpath
addpath trainingdata;
addpath processdata;
%% define path
currentPath = cd;

% input path
TrainingPath = fullfile(currentPath,'trainingdata') ;
TargetPath = fullfile(currentPath,'processdata') ; % path to files containing initial feature coordinates

%userName = 'tree80';
%userName = 'sciencebridge80';
userName = 'canoe60by80from800to1000';


%% Get training images
[fileNames, numImages] = gzget_training_images( TrainingPath, userName) ;
testImage = imread(fileNames{1});
if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

data = zeros([h*w, numImages], 'double');
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  if isrgb(currentImage)
    currentImage = rgb2gray(currentImage);
  end
  im = reshape(currentImage,w*h,1);
  data(:, fileIndex) = im(:,1);
end
DataTrain=data;

%% Get images to be processed
[fileNames2, numImages2] = gzget_training_images( TargetPath, userName) ;
testImage2 = imread(fileNames2{1});
if isrgb(testImage2)
    testImage2 = rgb2gray(testImage2);
end
[h, w] = size(testImage2);

data2 = zeros([h*w, numImages2], 'double');
for fileIndex2 = 1:numImages2
  currentImage2 = imread(fileNames2{fileIndex2});
  if isrgb(currentImage2)
    currentImage2 = rgb2gray(currentImage2);
  end
  im2 = reshape(currentImage2,w*h,1);
  data2(:, fileIndex2) = im2(:,1);
end

M=data2;
imSize=[h,w];
save(userName,'imSize','DataTrain','M');

