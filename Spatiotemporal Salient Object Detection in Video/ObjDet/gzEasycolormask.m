% Gao Zhi generate the color mask!
clc ; clear all; close all ;
% addpath

%% read first images
% testImage1 = imread(fileNamesGT{1});
testImage1 = imread('gt000981.png');
if isrgb(testImage1)
    testImage1 = rgb2gray(testImage1);
end
[h1, w1] = size(testImage1);

% testImage2 = imread(fileNamesOUR{1});
testImage2 = imread('bin000981.png');
if isrgb(testImage2)
    testImage2 = rgb2gray(testImage2);
end
[h2, w2] = size(testImage2);

%%
if ((h1 ~= h2)|| (w1 ~= w2))
    disp('Two sequence do not match! something error in the folder setting!') ;
end

%%
% for fileIndex = 1:numImagesGT
    dataR = zeros([h1, w1], 'double');
    dataG = zeros([h1, w1], 'double');
    dataB = zeros([h1, w1], 'double');
    dataRGB = zeros([h1,w1, 3], 'double');
    
%     GTimg = imread(fileNamesGT{fileIndex});
    GTimg = imread('gt000981.png');
    if isrgb(GTimg)
    GTimg = rgb2gray(GTimg);
    end
    BWGT = im2bw(GTimg,0.01);
    
    %GTimg = imread(fileNamesGT{fileIndex});
    %BWGT = im2bw(GTimg,0.1);
    gtmp1= GTimg>40;
    gtmp2= GTimg<40;
    gtmpunknow = gtmp1 & gtmp2;
    gtmpunknow = 1 - gtmpunknow;

    %OURmask = imread(fileNamesOUR{fileIndex});
    OURmask = imread('bin000981.png');
    if isrgb(OURmask)
    OURmask = rgb2gray(OURmask);
    end
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

    fileIndex=981;
    partnames = sprintf('zzColor%06d.bmp',fileIndex);
%     outputFileName  = strcat(userName, partnames);
%     outputFileNames = fullfile(destDir, outputFileName);
    imwrite(dataRGB, partnames);
% end