
% clc;
% clear all;
% close all;

% addpath
addpath(genpath('data'));
addpath(genpath('results'));

% methodchoose =1; % original RPCA code
methodchoose = 2; % block RPCA method

%% define images path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data');     % data\JHMDB\images;
flowPath = fullfile(currentPath,'flowdata');  % path to files containing initial flow info

% userName = 'April_09_brush_hair_u_nm_np1_ba_goo_0';
% userName = 'Aussie_Brunette_Brushing_Long_Hair_brush_hair_u_nm_np1_ba_med_3';
% userName = 'smtreeman100320';
% userName = 'raincar100';
userName = 'carheavyrain';
% userName = 'cardeemlight';
% userName = 'waterSurface';
% userName = 'walkStraight';
% userName = 'scienceblue';
% userName = 'carlodge1';
% userName = 'sciencenight';
% userName = 'scienceBigRain';
% userName = 'basicpart';
% userName = 'lightswitchpart';
% userName = 'noisynightpart';
% userName = 'whiteCar';
% userName = 'lightswithprehalf';
% userName = 'lightswitchselect';
% userName = 'lightswitchhalfhalf';
% userName = 'switchlightlab';
% userName = 'lightswitchshanmo';
% userName = 'boats';
% userName = 'flockpart';
% userName = 'flockselect';
% userName = 'manyBirds';
% userName = 'oneBirdFly';
% userName = 'peds';

% output path
if methodchoose==2
    % destRoot = fullfile(currentPath,'gzblock2and1adaptivethreshresults') ;
    destRoot = fullfile(currentPath,'gzblock2and1 1st round');
else  
    destRoot = fullfile(currentPath,'gzblockrpcaresults');
end

destDir = fullfile(destRoot,userName);

if ~exist(destDir,'dir')
    mkdir(destRoot,userName);
end

%% Get training images
[fileNames, numImages] = gzget_training_images(imagePath, userName);

%% read first image to get the size info
testImage = imread(fileNames{1});
if size(testImage, 3) > 1
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

%% read image file, to establish the data matrix
data = zeros([h*w, numImages], 'double');
for fileIndex = 1:numImages
    currentImage = imread(fileNames{fileIndex});
    if size(currentImage, 3) > 1
        currentImage = rgb2gray(currentImage);
    end
    im = reshape(currentImage,w*h,1);
    data(:, fileIndex) = im(:,1);
end
data=data/255.0;

%% read flow data, u v and norm %%
uflowdata = zeros([h*w, numImages-1], 'double');
vflowdata = zeros([h*w, numImages-1], 'double');
normflowdata = zeros([h*w, numImages-1], 'double');

flowDir= fullfile(flowPath, userName);
 
for j = 1:numImages - 1,
    infName0 = sprintf('Flowu%d.mat',j);
    infName00 = fullfile(flowDir, infName0);
    load(infName00,'u');
    utemp = reshape(u,w*h,1);
    uflowdata(:, j) = utemp(:,1);
    
    infName1 = sprintf('Flowv%d.mat',j);
    infName11 = fullfile(flowDir, infName1);
    load(infName11,'v');
    vtemp = reshape(v,w*h,1);
    vflowdata(:, j) = vtemp(:,1);
        
    infName2 = sprintf('norm%d.mat',j);
    %infName2 = sprintf('norm%d.mat',j);
    infName22 = fullfile(flowDir, infName2);
    load(infName22,'norm');
    normtemp = reshape(norm,w*h,1);
    normflowdata(:, j) = normtemp(:,1);
end 

%% the flow seems too noisy, norm is not stable, below functions are ok but useless
% gzFlowNormalizationShow(normflowdata);
% FlowPattern = gzFlowMagpattern(normflowdata, 4);

%% to run all kinds of RPCA
if methodchoose==2
% test save temporal data and load it
%     rtemp =rand(10,10);
%     save('blockinfo.mat', 'rtemp');
%     load('blockinfo.mat', 'rtemp');
% test how to save file in certain folder
%     gzTestSavefloder(data);

    % First-pass RPCA
    % A -> Background; E ->Moving objects (RPCA method)
    [A_hat, E_hat, iter] = gzblock2and1norm_inexact_alm_rpca(data);
    
    % value in MaskBlock1Pixle2 is 1,0.5,and 0. 1 is the block&pixel,0.5 is block mask    
    % value in E_blockmask is 1 or 0. represent block is maksked
    
    % for treeman case use 0.3; for raincar case use 0.2; for car case use 0.1; 
    
    % 1) blocksize = 8 case
    % [MaskBlock1Pixle2 E_blockmask]  = gzAnalyzeEmatrixFindMasks(E_hat,0.3);
    % PixelEasyMask  = gzEasyPixelMask(E_hat, 0.2);
    
    % 2) blocksize = 4 case
%     [MaskBlock1Pixle2 E_blockmask PixelMask]  = gzEmatrixFindMasksChangingblocksize(E_hat, 0.3, 4);
%     PixelEasyMask  = gzEasyPixelMask(E_hat, 0.2);
    
%     [MaskBlock1Pixle2 E_blockmask PixelMask]  = gzEmatrixFindMasksChangingblocksize(E_hat, 0.1, 4);
%     PixelEasyMask  = gzEasyPixelMask(E_hat, 0.05);
    
    % for carbluelight case
    [MaskBlock1Pixle2, E_blockmask, PixelMask] = gzEmatrixFindMasksChangingblocksize(E_hat, 0.1, 4);
    PixelEasyMask  = gzEasyPixelMask(E_hat, 0.11);
    
%     PixelMaskRelativeChange = gzComputeRelativeChangePixelmask(data, A_hat, 0.1);
%     [MaskBlock1Pixle2 E_blockmask] = gzEmatrixFindMasksChangingblocksize(PixelMaskRelativeChange, 0.2, 4);
    
    % Motion saliency estimation
    PixelSeedMask = gzMotionConsistentAnalysis(MaskBlock1Pixle2, E_blockmask, PixelMask, uflowdata, vflowdata, normflowdata, data);
    
    [TargetCoord, TargetMask] = gzGetTarget(PixelEasyMask,PixelSeedMask);    
%     save('blockinfo.mat', 'TargetCoord');
    gzShowTargetinImg(data, TargetCoord);
   
    % % Second-pass RPCA with the location and size of the likely outlier blocks estimated (E_hat)
    % clear A_hat E_hat MaskBlock1Pixle2 E_blockmask PixelMask PixelEasyMask PixelSeedMask;
    % [A_hat, E_hat, iter] = gzBlockRPCAwithTargetBlockinfo(data);
    
    [Mask_matrix, E_blockmask] = gzAnalyzeEmatrix(E_hat);
    [A_hat, E_hat, iter] = gzAdaptiveThreshblock2and1norm_inexact_alm_rpca(data, E_blockmask);
    PixelEasyMask = gzEasyPixelMask(E_hat, 0.25);
    [TargetCoord2, TargetPix]= gzGetTarget2(PixelEasyMask, TargetMask);
    gzShowTargetinImg(data, TargetCoord2);

end 

AANames = fullfile(destDir, 'Amatrix.mat'); 
save(AANames,'A_hat') ;
EENames = fullfile(destDir, 'Ematrix.mat'); 
save(EENames,'E_hat') ;

%% save the data as images (low rank and sparse matrix)
for fileIndex = 1:numImages
%     immask = Mask_matrix(:,fileIndex);
%     immask = reshape(immask,h,w);
%     immask = immask * 255;

%     outputMask  = sprintf('mask%d.jpg',fileIndex);
%     outputMasks = fullfile(destDir, outputMask); 
%     imwrite(immask, outputMasks);
    
    im1 = E_hat(:,fileIndex);
    for j = 1:h*w,
%         if (abs(im1(j,1))) < 0.02,
% 	        im1(j,1) = 0;
%         else im1(j,1)=255;
%         end
        
        if im1(j,1) < 0
            im1(j,1) = -1*im1(j,1);   
            % im1(j,1) = 1; % to much negative
        end
    end
    im1 = reshape(im1,h,w);
    % im1 = im1 * 255;
    % not display now
    % figure; imshow(uint8(im1));
    
    outputFileNameLowrank  = sprintf('Etarget%d.jpg',fileIndex);
    outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
    imwrite(im1, outputFileNamesLowrank);
    % imwrite(uint8(im1), outputFileNamesLowrank);  % if *255
    
    im2 = A_hat(:,fileIndex);
    im2 = reshape(im2,h,w);
    im2 = im2 * 255;
    % not display now
    % figure; imshow(uint8(im2));
    
    outputFileNameSparse  = sprintf('Abackground%d.jpg',fileIndex);
    outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
    % imwrite(im2, outputFileNamesSparse);
    imwrite(uint8(im2), outputFileNamesSparse);
end
