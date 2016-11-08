% Gao Zhi generate the color mask!
clc ; clear all; close all ;
% addpath
addpath groundtruthmask;
addpath ourresultmask;
addpath colormask;

%% define images' path
currentPath = cd;
% input path
imagePathgroundtruth = fullfile(currentPath,'groundtruthmask') ;
imagePathourresult   = fullfile(currentPath,'ourresultmask') ;

%userName = 'aligned00900aligned295after00120to00210';
%userName = 'aligned00900aligned001to294range00085to00157';
% userName = 'sidewalkallspecialpart2';
%userName = 'original00900onwardsaligned295afterRASL';
%userName = 'original00900onwardsaligned001to294RASL';
% userName = 'sidewalkallspecialpart2RASL';
% userName = 'aligned00900aligned295after00120to00210DECOLOR';
%userName = 'park240by440to540MADDALENA';
%userName = 'park240by440to540BARNICH';
%userName = 'park240by440to540KDE';
% userName = 'park240by440to540GMM';
% userName = 'park120by440to540DECOLOR';
% userName = 'park240by440to540pcp';
userName = 'park240by440to540our';







% output path
destRoot = fullfile(currentPath,'colormask') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNamesGT, numImagesGT] = gzget_training_images( imagePathgroundtruth, userName) ;
[fileNamesOUR, numImagesOUR] = gzget_training_images( imagePathourresult, userName) ;

%% read first images
testImage1 = imread(fileNamesGT{1});
if isrgb(testImage1)
    testImage1 = rgb2gray(testImage1);
end
[h1, w1] = size(testImage1);

testImage2 = imread(fileNamesOUR{1});
if isrgb(testImage2)
    testImage2 = rgb2gray(testImage2);
end
[h2, w2] = size(testImage2);

%%
if ((numImagesGT ~= numImagesOUR) || (h1 ~= h2)|| (w1 ~= w2))
    disp('Two sequence do not match! something error in the folder setting!') ;
end

%%
for fileIndex = 1:numImagesGT
    dataR = zeros([h1, w1], 'double');
    dataG = zeros([h1, w1], 'double');
    dataB = zeros([h1, w1], 'double');
    dataRGB = zeros([h1,w1, 3], 'double');
    
    GTimg = imread(fileNamesGT{fileIndex});
    BWGT = im2bw(GTimg,0.1);
    
    %GTimg = imread(fileNamesGT{fileIndex});
    %BWGT = im2bw(GTimg,0.1);
    gtmp1= GTimg>160;
    gtmp2= GTimg<200;
    gtmpunknow = gtmp1 & gtmp2;
    gtmpunknow = 1 - gtmpunknow;

    OURmask = imread(fileNamesOUR{fileIndex});
    BWOUR = im2bw(OURmask,0.1);

    correctMask = BWGT & BWOUR;
    dataR       =correctMask*1.0;
    dataG       =correctMask*1.0;
    dataB       =correctMask*1.0;

    missMask    = BWGT & (~BWOUR);
    dataR       = dataR + missMask*1.0;

    falsealarmMask = (~BWGT) & BWOUR;
    dataB       = dataB + falsealarmMask*1.0;
    
    dataRGB(:,:, 1) = dataR.*gtmpunknow;
    dataRGB(:,:, 2) = dataG.*gtmpunknow;
    dataRGB(:,:, 3) = dataB.*gtmpunknow;

    partnames = sprintf('Color%06d.bmp',fileIndex);
    outputFileName  = strcat(userName, partnames);
    outputFileNames = fullfile(destDir, outputFileName);
    imwrite(dataRGB, outputFileNames);
end