%when there is too much training data, try to delete some usefulness

% clear
clc ; clear all; close all ;

addpath dataoftraining ;
addpath dataROImask ;
addpath dataoftrainingROI;


currentPath = cd;
% input path
imagePath = fullfile(currentPath,'dataoftraining') ;
%userName  = 'cam2new';
%userName  = 'cam4new';
userName  = 'cam1new';
maskNmae = 'generatecam1ROI.PNG';

maskPath = fullfile(currentPath,'dataROImask') ;
maskfile = fullfile(maskPath,maskNmae) ;
maskImage= imread(maskfile);
maskBinary=maskImage>250;
%imshow(maskBinary);


destRoot = fullfile(currentPath,'dataoftrainingROI');
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

for fileIndex = 1 : numImages
    currentImage = imread(fileNames{fileIndex});
    
    rstImg(:,:,1) = double(currentImage(:,:,1)).*maskBinary;
    rstImg(:,:,2) = double(currentImage(:,:,2)).*maskBinary;
    rstImg(:,:,3) = double(currentImage(:,:,3)).*maskBinary;
    rstImg        =rstImg/255.0;
    %imshow(rstImg);
    
    outputname  = sprintf('trainingROI%d.jpg',fileIndex);
    outputnames = fullfile(destDir, outputname);
    imwrite(rstImg, outputnames);
    %imwrite(uint8(BkImage), outputFileNamesSparse);
end 
  