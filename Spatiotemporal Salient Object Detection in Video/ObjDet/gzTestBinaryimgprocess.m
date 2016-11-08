clc ;
clear all;
close all ;

% addpath
addpath maskori ;
addpath maskrefine ;

currentPath = cd;
% input path
imagePath = fullfile(currentPath,'maskori') ;
userName = 'Frontier_Cam01testaccumulate';
% userName = 'Frontier_Cam01testaccumulate38to73';


%
destRoot = fullfile(currentPath,'maskrefine') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end


[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

for fileIndex = 1:numImages
    %     currentImage = imread(fileNames{fileIndex});
    %     IMin0 = imread('gzBinaryMask000001.bmp');
    
    IMin0 = imread(fileNames{fileIndex});
    % Masktmp = bwmorph(IMin0, 'clean');
    Masktmp = bwmorph(IMin0, 'erode');
    Masktmp = bwmorph(Masktmp, 'dilate');
    % Masktmp = bwmorph(Masktmp, 'fill');
    Masktmp = bwmorph(Masktmp, 'dilate');
    %     imshow(IMin0);
    %
    %     figure, imshow(Masktmp);
    %     % dataR0=IMin0(:,:,1);

    outputFileName = sprintf('maskrefine%06d.bmp',fileIndex);
    outputFileNames = fullfile(destDir, outputFileName);
    imwrite(Masktmp, outputFileNames);
end
