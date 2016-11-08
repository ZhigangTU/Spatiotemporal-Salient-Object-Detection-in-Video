
% addpath
addpath(genpath('data'));
addpath(genpath('results'));
addpath PROPACK;

% methodchoose =1; % original RPCA code
methodchoose = 2;  % block RPCA method


%% define images path
currentPath = cd;

% input path
userName = '012';
% imagePath = fullfile(currentPath,'data\input\Diving-Side');         % JHMDB\images\climb_stairs; Diving-Side
% flowPath = fullfile(currentPath,'flowdata\Diving-Side');            % path to files containing initial flow info

imagePath = fullfile('E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images','SkateBoarding-Front');  % JHMDB\images\climb_stairs; data\input\Diving-Side
flowPath = fullfile('E:\Action Recognition\Pose-CNN\cache\OFUCFSport','SkateBoarding-Front');         % flowdata\Diving-Side

outPath = fullfile(currentPath,'data\output\SkateBoarding-Front',userName); % JHMDB\brush_hair

if ~exist(outPath,'dir')
    mkdir(outPath); 
end

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
framesdata = cell(numImages, 1);
% framesName = cell(numImages, 1);
for fileIndex = 1:numImages
    currentImage = imread(fileNames{fileIndex});
    framesdata{fileIndex} = double(currentImage);
%     [~, frameName] = fileparts(fileNames{fileIndex});
%     framesName{fileIndex} = frameName;
    if size(currentImage, 3) > 1
        currentImage = rgb2gray(currentImage);
    end
    im = reshape(currentImage,w*h,1);
    data(:, fileIndex) = im(:,1);
end
data = data/255.0;

%% read flow data, u v and norm %%
uflowdata = zeros([h*w, numImages-1], 'double');
vflowdata = zeros([h*w, numImages-1], 'double');
flowdata = cell(1, numImages-1);
magnitude = cell(1, numImages-1);
magflowdata = zeros([h*w, numImages-1], 'double');

flowDir = fullfile(flowPath, userName);
 
for j = 1:numImages - 1,
    infName0 = sprintf('flowu%d.mat',j);
    infName00 = fullfile(flowDir, infName0);
    load(infName00,'u');
    uflow = u;
    utemp = reshape(u,w*h,1);
    uflowdata(:, j) = utemp(:,1);
    
    infName1 = sprintf('flowv%d.mat',j);
    infName11 = fullfile(flowDir, infName1);
    load(infName11,'v');
    vflow = v;
    vtemp = reshape(v,w*h,1);
    vflowdata(:, j) = vtemp(:,1);
        
%     infName2 = sprintf('norm%d.mat',j);
%     % infName2 = sprintf('norm%d.mat',j);
%     infName22 = fullfile(flowDir, infName2);
%     load(infName22,'norm');
%     magtemp = reshape(norm,w*h,1);
%     magflowdata(:, j) = magtemp(:,1);
    
    flowdata{j} = cat(3,uflow,vflow);
    Flowgradient = getFlowGradient(flowdata{j});
    magnitude{j} = getMagnitude(Flowgradient);
    magflowdata(:, j) = reshape(magnitude{j},w*h,1);
end 


% %% Performing our Appearance-Motion based object detection
% % Video Saliency Detection (Appearance + Motion)-->Probability + InOutMaps
% [frameIOMap,frameObject,frameObjectBoxComb, frameEnergy] = SaliencyObject(flowdata,magnitude,framesdata,h,w,numImages,outPath);  % frameEnergy-->superpixels+frameObject
% % TShowTargetinImg(frameEnergy, outPath, h, w);


%% Performing the RPCA
if methodchoose==2
% test save temporal data and load it
%     rtemp =rand(10,10);
%     save('blockinfo.mat', 'rtemp');
%     load('blockinfo.mat', 'rtemp');
% test how to save file in certain folder
%     gzTestSavefloder(data);

    % 1. First-pass RPCA
    % A -> Background; E ->Moving objects (RPCA method)
%     [A_hat, E_hat, iter] = gzblock2and1norm_inexact_alm_rpca(data);  
    [A_hat, E_hat, iter] = Block2and1norm_inexact_alm_rpca(data, h, w);
    
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
%     [MaskBlock1Pixle2, E_blockmask, PixelMask] = gzEmatrixFindMasksChangingblocksize(E_hat, 0.1, 4);
%     PixelEasyMask  = gzEasyPixelMask(E_hat, 0.11);
    [MaskBlock1Pixle2, E_blockmask, PixelMask, E_blockmaskSize] = TEmatrixFindMasksChangingblocksize(E_hat, 0.05, 4, h, w); % nSetBlock=[4, 8]
    % 0.05; 0.1; 0.11; 0.2
    PixelEasyMask = TEasyPixelMask(E_hat, 0.05, h, w);  % Large(0.2)-->more edge features (slow), Small(0.09)-->few edge features (fast)
       
%     PixelMaskRelativeChange = gzComputeRelativeChangePixelmask(data, A_hat, 0.1);
%     [MaskBlock1Pixle2 E_blockmask] = gzEmatrixFindMasksChangingblocksize(PixelMaskRelativeChange, 0.2, 4);
    
    % 2. Motion saliency estimation
%     PixelSeedMask = gzMotionConsistentAnalysis(MaskBlock1Pixle2, E_blockmask, PixelMask, uflowdata, vflowdata, normflowdata, data);
    PixelSeedMask = TMotionConsistentAnalysis(MaskBlock1Pixle2, E_blockmask, PixelMask, uflowdata, vflowdata, magflowdata, data, h, w, E_blockmaskSize);
    
%     [TargetCoord, TargetMask] = TGetTarget(PixelEasyMask, PixelSeedMask, h, w, outPath);  % TargetMask --> logical labels  
    frameTargetBox = TargetBoxLab(framesdata, PixelEasyMask, PixelSeedMask, h, w, outPath);
    
    TargetObjects = TargetRPCAMotionFuse(frameTargetBox, frameObjectBoxComb, numImages, framesdata, magnitude, outPath);
    
    TargetObjectsCoord = sprintf('%s/TargetObjectsCoord.mat',outPath);
    save(TargetObjectsCoord,'TargetObjects');

    
%     % % 3. Second-pass RPCA with the location and size of the likely outlier blocks estimated (E_hat)
%     % clear A_hat E_hat MaskBlock1Pixle2 E_blockmask PixelMask PixelEasyMask PixelSeedMask;
%     % [A_hat, E_hat, iter] = gzBlockRPCAwithTargetBlockinfo(data);
%     
%     [Mask_matrix, E_blockmask] = TAnalyzeEmatrix(E_hat, h, w); % nSetBlock=[4, 8]
%     [A_hat, E_hat, iter] = TAdaptiveThreshblock2and1norm_inexact_alm_rpca(data, E_blockmask, h, w); % nSetBlock=[4, 8]
%     % 0.05; 0.1; 0.11; 0.20; 0.25
%     PixelEasyMask = TEasyPixelMask(E_hat, 0.25, h, w);
%     Time = 2; % second time
%     [TargetCoord2, TargetPix]= TGetTarget(PixelEasyMask, PixelSeedMask, h, w, outPath, Time); % TargetPix --> logical labels after Second-pass RPCA
%     gzShowTargetinImg(data, TargetCoord2, userName, h, w, Time);
    
end 
