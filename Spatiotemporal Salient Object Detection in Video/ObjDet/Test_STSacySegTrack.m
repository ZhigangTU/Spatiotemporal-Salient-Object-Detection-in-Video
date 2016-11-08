function Test_STSacySegTrack(videoName)
% videoName = 'Riding-Horse';
% SeqName   = '001';

addpath(genpath('data'));
addpath(genpath('results'));
addpath PROPACK;

methodchoose = 2;  % 2-->block RPCA method; % 1-->original RPCA code

%% define images path
currentPath = cd;

%% input path
imagePath = fullfile('F:\SegTrackV1\Images');  % JHMDB\images\climb_stairs; data\input\Diving-Side
flowPath  = fullfile('E:\Action Recognition\Pose-CNN\cache\OFSegTrackV1-Brox2011', videoName);        % flowdata\Diving-Side

outPath = fullfile(currentPath,'data\output\SegTrackV1',videoName);   % JHMDB\brush_hair
outPathRPCA = fullfile(currentPath,'data\output\SegTrackV1',videoName,'RPCA');  
outPathFOS  = fullfile(currentPath,'data\output\SegTrackV1',videoName,'FOS');   

if ~exist(outPath,'dir')
    mkdir(outPath); 
end
if ~exist(outPathRPCA,'dir')
    mkdir(outPathRPCA); 
end
if ~exist(outPathFOS,'dir')
    mkdir(outPathFOS); 
end

% output path
if methodchoose==2
    % destRoot = fullfile(currentPath,'gzblock2and1adaptivethreshresults') ;
    destRoot = fullfile(currentPath,'gzblock2and1 1st round');
else  
    destRoot = fullfile(currentPath,'gzblockrpcaresults');
end

destDir = fullfile(destRoot);

if ~exist(destDir,'dir')
    mkdir(destRoot);
end

%% Get training images
[fileNames, numImages] = gzget_training_images(imagePath, videoName);

%% read first image to get the size info
testImage = imread(fileNames{1});
if size(testImage, 3) > 1
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

%% read image file, to establish the data matrix
data = zeros([h*w, numImages], 'double');
framesdata = cell(numImages, 1);  %framesName = cell(numImages, 1);
for fileIndex = 1:numImages
    currentImage = imread(fileNames{fileIndex});
    framesdata{fileIndex} = currentImage;
    if size(currentImage, 3) > 1
        currentImage = rgb2gray(currentImage);
    end
    im = reshape(currentImage,w*h,1);
    data(:, fileIndex) = im(:,1);
end
data = data/255.0;

%% read flow data, u v and norm %%
uOrgflowdata = zeros([h*w, numImages-1], 'double');
vOrgflowdata = zeros([h*w, numImages-1], 'double');
magflowdata = zeros([h*w, numImages-1], 'double');
flowdata    = cell(numImages-1,1);
flowOrgdata = cell(numImages-1,1);
magnitude   = cell(numImages-1,1);
flowDir = flowPath;
shot = 1;

flowmethod = 'broxPAMI2011';
flowdata = loadFlow(flowDir, shot, flowmethod);
flowOrgdata = [];
for j = 1:numImages - 1
    Flowgradient = getFlowGradient(flowdata{j});
    magnitude{j} = getMagnitude(Flowgradient);
    magflowdata(:, j) = reshape(magnitude{j},w*h,1);
end

% Load GroundTruth
options.infolder = fullfile(imagePath, videoName);
options.SeqName = videoName;
gtFile = fullfile(options.infolder, 'GroundTruth', sprintf('GT%s.mat', videoName)); %'1'
groundTruth = load( gtFile );
groundTruth = groundTruth.groundTruth; 

% for j = 1:numImages - 1
%     % 2) Load Brox or classic-NL
%     infName0 = sprintf('flowu%d.mat',j);
%     infName00 = fullfile(flowDir, infName0);
%     load(infName00,'u');    
%     utemp = reshape(u,w*h,1);
%     uOrgflowdata(:, j) = utemp(:,1);
%     uIntflow = int16(u);
%     
%     infName1 = sprintf('flowv%d.mat',j);
%     infName11 = fullfile(flowDir, infName1);
%     load(infName11,'v');    
%     vtemp = reshape(v,w*h,1);
%     vOrgflowdata(:, j) = vtemp(:,1);
%     vIntflow = int16(v);
%     
%     flowOrgdata{j} = cat(3,u,v);
% 
%     flowdata{j}  = cat(3,uIntflow,vIntflow);
%     Flowgradient = getFlowGradient(flowdata{j});
%     magnitude{j} = getMagnitude(Flowgradient);
%     magflowdata(:, j) = reshape(magnitude{j},w*h,1);
    
%     % 3) Load MDP flow
%     for j = 1:nframe - 1
%         flowName = sprintf('flow%04d_mdpof.flo',j);
%         flowPath = fullfile(flowDir, flowName);
%         uv = readFlowFile(flowPath);
% 
%         flowOrgdata{j} = cat(3,uv(:,:,1),uv(:,:,2));
% 
%         uflow = int16(uv(:,:,1));
%         vflow = int16(uv(:,:,2));
%         flowdata{j}  = cat(3,uflow,vflow);
%         Flowgradient = getFlowGradient(flowOrgdata{j});
%         magnitude{j} = getMagnitude(Flowgradient);
%         magflowdata(:, j) = reshape(magnitude{j},w*h,1);
%     end 
% end 


%% Performing our Appearance-Motion based object detection --> Video Saliency Detection (Appearance + Motion)

%% 1) Performing the FOS
addpath(genpath('FastVideoSegment'));
TData.Orgflow= flowOrgdata;
TData.flow   = flowdata;
TData.imgs   = framesdata;
TData.flowmag= magnitude;

segmfolder = fullfile(outPathFOS, 'segmentations', 'VideoRapidSegment');
if(~exist(segmfolder, 'dir'))
    mkdir(segmfolder)
end
FOSSegLabelsName = fullfile(segmfolder, 'segmentation.mat');
if exist(FOSSegLabelsName, 'file')
    FOSSegLabels = load(fullfile(segmfolder, 'segmentation.mat'));
    frameFOSSegLabels = FOSSegLabels.FOSsegmentation;
else 
    frameFOSSegLabels = FastObjSeg(TData,numImages,flowPath,imagePath,outPathFOS,videoName,segmfolder);
end
% %segmentationLarg
% FOSSegLabelsLargName = fullfile(segmfolder, 'segmentationLarg.mat');
% if exist(FOSSegLabelsLargName, 'file')
%     FOSSegLargLabels = load(FOSSegLabelsLargName);
%     frameFOSSegLabels = FOSSegLargLabels.segmentationLarg;
% else 
%     frameFOSSegLabels = FastObjSeg(TData,numImages,flowPath,imagePath,outPathFOS,videoName,segmfolder);
% end

%% 2) Performing the RPCA
% test save temporal data and load it
%     rtemp =rand(10,10);
%     save('blockinfo.mat', 'rtemp');
%     load('blockinfo.mat', 'rtemp');
% test how to save file in certain folder
%     gzTestSavefloder(data);
RPCATargetLabelsName = fullfile(outPathRPCA,'001-001-5-20','TargetObjectsLabels.mat');
if exist(RPCATargetLabelsName, 'file')
    TargetLabels = load(RPCATargetLabelsName);
    frameTargetLabels = TargetLabels.frameTargetLabels;
else
    % 1. First-pass RPCA
    % A -> Background; E ->Moving objects (RPCA method)
%     [A_hat, E_hat, iter] = gzblock2and1norm_inexact_alm_rpca(data);  
    [A_hat, E_hat, iter] = Block2and1norm_inexact_alm_rpca(data, h, w);

    % value in MaskBlock1Pixle2 is 1,0.5,and 0. 1 is the block&pixel,0.5 is block mask    
    % value in E_blockmask is 1 or 0. represent block is maksked

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
    [MaskBlock1Pixle2, E_blockmask, PixelMask, E_blockmaskSize] = TEmatrixFindMasksChangingblocksize(E_hat, 0.02, 4, h, w); % nSetBlock=[4, 8];0.1-->background
    % 0.05; 0.1; 0.11; 0.2
    PixelEasyMask = TEasyPixelMask(E_hat, 0.02, h, w);  % Large(0.2)-->more edge features (slow), Small(0.09)-->few edge features (fast);(UCF-->0.13;SegTrack-->0.10)

%     PixelMaskRelativeChange = gzComputeRelativeChangePixelmask(data, A_hat, 0.1);
%     [MaskBlock1Pixle2 E_blockmask] = gzEmatrixFindMasksChangingblocksize(PixelMaskRelativeChange, 0.2, 4);

    % 2. Motion saliency estimation
%     PixelSeedMask = gzMotionConsistentAnalysis(MaskBlock1Pixle2, E_blockmask, PixelMask, uOrgflowdata, vOrgflowdata, normflowdata, data);
    PixelSeedMask = TMotionConsistentAnalysis(MaskBlock1Pixle2, E_blockmask, PixelMask, uOrgflowdata, vOrgflowdata, magflowdata, data, h, w, E_blockmaskSize);

%     [TargetCoord, TargetMask] = TGetTarget(PixelEasyMask, PixelSeedMask, h, w, outPath);  % TargetMask --> logical labels  

    [frameTargetBox, frameTargetLabels] = TargetBoxLab(framesdata, PixelEasyMask, PixelSeedMask, h, w, outPathRPCA);  %TargetObjects=TargetRPCAMotionFuse(frameTargetBox,frameObjectBoxComb,numImages);
      
    % Save output
    TargetObjectsBoxs = sprintf('%s/TargetObjectsBoxs.mat',outPathRPCA);
    save(TargetObjectsBoxs,'frameTargetBox');
    TargetObjectsLabels = sprintf('%s/TargetObjectsLabels.mat',outPathRPCA);
    save(TargetObjectsLabels,'frameTargetLabels');
    
    RPCAsegmentationLarg = getLargestSegmentAndNeighbours(frameTargetLabels);
    for idx = 1:length(frameTargetLabels)
        frameRPCASeg = frameTargetLabels{idx};
        outObjectPath = sprintf('%s/frameRPCASegFram%03d.jpg',outPathRPCA, idx);
        imwrite(frameRPCASeg, outObjectPath);
        
        frameRPCASegLarg = RPCAsegmentationLarg{idx};
        outObjectPath = sprintf('%s/frameRPCASegLargFram%03d.jpg',outPathRPCA, idx);    
        imwrite(frameRPCASegLarg, outObjectPath);
    end
    TargetObjectsLargLabels = sprintf('%s/TargetObjectsLargLabels.mat',outPathRPCA);
    save(TargetObjectsLargLabels,'RPCAsegmentationLarg');

    avgMislabelled = getAverageMislabelledPixels(options, frameTargetLabels);
    fprintf( 'Average number of mislabelled pixels for %s: %i\n', videoName, avgMislabelled );
    avgMislabelled = getAverageMislabelledPixels(options, RPCAsegmentationLarg);
    fprintf( 'Average number of mislabelled pixels for %s: %i\n', videoName, avgMislabelled );
end
        
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
%     gzShowTargetinImg(data, TargetCoord2, SeqName, h, w, Time);


%%  Salient object detection via objectness measure? Foreground Connectivity
addpath('Foreground Connectivity\Funcs');

%% 1. Parameter Settings
doFrameRemoving = false; % true; false;
useSP = true;            % You can set useSP = false to use regular grid for speed consideration

%% 2. Saliency Map Calculation
% (1)RPCA
framewCtr = cell(numImages, 1);
frameoptwCtr = cell(numImages, 1);
frameoptFgConSacyMapsPCA = cell(numImages, 1);
frameRefinConSacyMapsPCA = cell(numImages, 1);
% (2)FOS
framewCtrFOS = cell(numImages, 1);
frameoptwCtrFOS = cell(numImages, 1);
frameoptFgConSacyMapsFOS = cell(numImages, 1);
frameRefinConSacyMapsFOS = cell(numImages, 1);
% (3) Combined Spatio-temperal Saliency Maps
frameSpaTempConSacyMaps = cell(numImages, 1);
% frameSpaTempConSacyMapsLarg = cell(numImages, 1);

framePixelList = cell(numImages, 1);
for j = 1:numImages
    % Pre-Processing: Remove Image Frames
    srcImg = uint8(framesdata{j});
    if doFrameRemoving
        [noFrameImg, frameRecord] = removeframe(srcImg, 'sobel');
        [h, w, chn] = size(noFrameImg);
    else
        noFrameImg = srcImg;
        [h, w, chn] = size(noFrameImg);
        frameRecord = [h, w, 1, h, 1, w];
    end
    
    %% Segment input rgb image into patches (SP/Grid)
    pixNumInSP = 150;  %pixels in each superpixel:[400,500,600,700,800]
    spnumber = round( h * w / pixNumInSP );     %super-pixel number for current image    
    if useSP
        [idxImg, adjcMatrix, pixelList] = SLIC_Split(noFrameImg, spnumber);
    else
        [idxImg, adjcMatrix, pixelList] = Grid_Split(noFrameImg, spnumber);        
    end
    framePixelList{j} = pixelList;
    
    %% Get super-pixel properties
    spNum = size(adjcMatrix, 1);
    meanRgbCol = GetMeanColor(noFrameImg, pixelList);
    meanLabCol = colorspace('Lab<-', double(meanRgbCol)/255);
    meanPos = GetNormedMeanPos(pixelList, h, w);
    bdIds = GetBndPatchIds(idxImg);           %super-pixels on image boundary
    colDistM = GetDistanceMatrix(meanLabCol); %pair-wise distance matrix between each rows in LabCol feature
    posDistM = GetDistanceMatrix(meanPos);    %pair-wise distance matrix between each rows in spatial
    [clipVal, geoSigma, neiSigma] = EstimateDynamicParas(adjcMatrix, colDistM);
    
    %% Saliency --> Foreground Connectivity
    [bgProb, bdCon, bgWeight] = EstimateBgProb(colDistM, adjcMatrix, bdIds, clipVal, geoSigma);
    % 1) Computing Foreground Weight (Wfg) according to RPCA
    [fgwCtrPCA, TwCtr] = SaliencyObjectnessTu(frameTargetLabels{j},h,w,idxImg,pixelList,adjcMatrix,colDistM,clipVal);  % srcName = sprintf('%05d',j); srcName(1:end-4)
    NonwCtrPCA = isnan(fgwCtrPCA);
    if sum(NonwCtrPCA(:))>1
        if j > 1
            framewCtr{j} = framewCtr{j-1};
            frameoptwCtr{j} = frameoptwCtr{j-1};
        else
            framewCtr{j} = zeros(size(fgwCtrPCA))+eps;
            frameoptwCtr{j} = zeros(size(fgwCtrPCA))+eps;
        end
        optFgConSacyMapsPCA = frameTargetLabels{j};
        frameoptFgConSacyMapsPCA{j} = frameTargetLabels{j};  
    else
        optwCtrPCA = SaliencyOptimization(adjcMatrix, bdIds, colDistM, neiSigma, bgWeight, fgwCtrPCA); %fgwCtr-->foreground weighted contrast 
        noSuffixName = sprintf('%05d',j);
        smapName = fullfile(outPath, strcat(noSuffixName, '_SObj.png'));
        optFgConSacyMapsPCA = SaveSaliencyMapTu(optwCtrPCA, pixelList, frameRecord, smapName, true);  

        framewCtr{j}    = fgwCtrPCA;
        frameoptwCtr{j} = optwCtrPCA;       
        frameoptFgConSacyMapsPCA{j} = optFgConSacyMapsPCA;
    end
    % Save the PCA SacyMapsPCA
    outObjectPath = sprintf('%s/FgConSacyMapsPCA%05d.jpg',outPath, j);    
    imwrite(optFgConSacyMapsPCA, outObjectPath);
    % Save the threshold Binary SacyMaps
    TMapsPCA = graythresh(optFgConSacyMapsPCA);
    outObjectPath = sprintf('%s/FgConSacyMapsPCASegTA%05d.jpg',outPath, j);
    imwrite(optFgConSacyMapsPCA>TMapsPCA, outObjectPath);
    
    % 2) Computing Foreground Weight (Wfg) according to Fast Object Segmentation on Optical flow (FOS)
    [fgwCtrFOS, TwCtrFOS] = SaliencyObjectnessTu(frameFOSSegLabels{j},h,w,idxImg,pixelList,adjcMatrix,colDistM,clipVal);
    NonwCtrFOS = isnan(fgwCtrFOS);    
    if sum(NonwCtrFOS(:))>1
        if j > 1
            framewCtrFOS{j} = framewCtrFOS{j-1};
            frameoptwCtrFOS{j} = frameoptwCtrFOS{j-1};  
        else
            framewCtrFOS{j} = zeros(size(fgwCtrFOS))+eps;
            frameoptwCtrFOS{j} = zeros(size(fgwCtrFOS))+eps;
        end
        optFgConSacyMapsFOS = frameFOSSegLabels{j};
        frameoptFgConSacyMapsFOS{j} = frameFOSSegLabels{j};   
    else
        optwCtrFOS = SaliencyOptimization(adjcMatrix, bdIds, colDistM, neiSigma, bgWeight, fgwCtrFOS);
        noSuffixName = sprintf('%05d',j);
        smapName = fullfile(outPath, strcat(noSuffixName, '_FOSObj.png'));
        optFgConSacyMapsFOS = SaveSaliencyMapTu(optwCtrFOS, pixelList, frameRecord, smapName, true);  

        framewCtrFOS{j}    = fgwCtrFOS;
        frameoptwCtrFOS{j} = optwCtrFOS;       
        frameoptFgConSacyMapsFOS{j} = optFgConSacyMapsFOS;
    end
    % Save the FOS SacyMapsFOS
    outObjectPath = sprintf('%s/FgConSacyMapsFOS%05d.jpg', outPath, j);    
    imwrite(optFgConSacyMapsFOS, outObjectPath);
    % Save the threshold Binary SacyMaps
    TMapsFOS = graythresh(optFgConSacyMapsFOS);
    outObjectPath = sprintf('%s/FgConSacyMapsFOSSegTA%05d.jpg',outPath, j);
    imwrite(optFgConSacyMapsFOS>TMapsFOS, outObjectPath);
    
    
    %% Spatio-temperal Fusion Refinement (Appearance-Motion)
    [SpaTempConSacyMap,RefinConSacyMapsPCA,RefinConSacyMapsFOS] = SacyMapsFuse(optFgConSacyMapsPCA,optFgConSacyMapsFOS,TMapsPCA,TMapsFOS);
    
    frameRefinConSacyMapsPCA{j} = RefinConSacyMapsPCA;    
    frameRefinConSacyMapsFOS{j} = RefinConSacyMapsFOS;
    frameSpaTempConSacyMaps{j}  = SpaTempConSacyMap;
 
    outObjectPath = sprintf('%s/SpaTempSacyMaps%05d.jpg', outPath, j); 
    imwrite(SpaTempConSacyMap, outObjectPath);    
    TSpaTmpMaps = graythresh(SpaTempConSacyMap);
    outObjectPath = sprintf('%s/SpaTempSacyMapsSegTAu%05d.jpg', outPath, j); 
    imwrite(SpaTempConSacyMap>TSpaTmpMaps, outObjectPath);
    outObjectPath = sprintf('%s/SpaTempSacyMapsSegTAu06%05d.jpg', outPath, j); 
    imwrite(SpaTempConSacyMap>TSpaTmpMaps*0.6, outObjectPath);
    
    outObjectPath = sprintf('%s/RefConSacyMapsPCA%05d.jpg', outPath, j); 
    imwrite(RefinConSacyMapsPCA, outObjectPath);
    outObjectPath = sprintf('%s/RefConSacyMapsFOS%05d.jpg', outPath, j); 
    imwrite(RefinConSacyMapsFOS, outObjectPath);
    
    if j == numImages
        figure; imshow(SpaTempConSacyMap)
    end
end
% Save Results
ConSacyMapsPCA = sprintf('%s/frameConSacyMapsPCA.mat',outPath);
save(ConSacyMapsPCA,'frameoptFgConSacyMapsPCA');
ConSacyMapsFOS = sprintf('%s/frameConSacyMapsFOS.mat',outPath);
save(ConSacyMapsFOS,'frameoptFgConSacyMapsFOS');
% Evaluate SacyMapsPCA
[prec, rec, mae] = evaluate_PR_MAE(frameoptFgConSacyMapsPCA, groundTruth);
fprintf('\noptFgConSacyMapsPCA: Precision=%f Recall=%f MAE=%f\n', prec, rec, mae);  
% Evaluate SacyMapsFOS
[prec, rec, mae] = evaluate_PR_MAE(frameoptFgConSacyMapsFOS, groundTruth);
fprintf('\noptFgConSacyMapsFOS: Precision=%f Recall=%f MAE=%f\n', prec, rec, mae); 
    
ReConSacyMapsPCA = sprintf('%s/frameRefConSacyMapsPCA.mat',outPath);
save(ReConSacyMapsPCA,'frameRefinConSacyMapsPCA');
ReConSacyMapsFOS = sprintf('%s/frameRefConSacyMapsFOS.mat',outPath);
save(ReConSacyMapsFOS,'frameRefinConSacyMapsFOS');
% Evaluate RefinSacyMapsPCA
[prec, rec, mae] = evaluate_PR_MAE(frameRefinConSacyMapsPCA, groundTruth);
fprintf('\nRefinoptFgConSacyMapsPCA: Precision=%f Recall=%f MAE=%f\n', prec, rec, mae);  
% Evaluate RefinSacyMapsFOS
[prec, rec, mae] = evaluate_PR_MAE(frameRefinConSacyMapsFOS, groundTruth);
fprintf('\nRefinoptFgConSacyMapsFOS: Precision=%f Recall=%f MAE=%f\n', prec, rec, mae); 

SpaTempSacyMaps = sprintf('%s/frameSpaTempSacyMaps.mat',outPath);
save(SpaTempSacyMaps,'frameSpaTempConSacyMaps');
% Evaluate fused SpaTempConSacyMap
[prec, rec, mae] = evaluate_PR_MAE(frameSpaTempConSacyMaps, groundTruth);
fprintf('\nFused-SpaTempConSacyMap: Precision=%f Recall=%f MAE=%f\n', prec, rec, mae);

SpixelList = sprintf('%s/SpixelList.mat',outPath);
save(SpixelList,'framePixelList');


%% Spatio-temperal Optimization
% % Load frameSpaTempConSacyMaps if it exist
% STSacyMapsLd = 'SpaTempSacy-002-002-5-20';
% frameSpaTempSacyMapsSav = fullfile(outPath, STSacyMapsLd, 'frameSpaTempSacyMaps.mat');
% if exist(frameSpaTempSacyMapsSav, 'file')    
%     EstSpaTempSacyMaps = load(frameSpaTempSacyMapsSav);
%     frameSpaTempConSacyMaps = EstSpaTempSacyMaps.frameSpaTempConSacyMaps;
% end
% [prec, rec, mae] = evaluate_PR_MAE(frameSpaTempConSacyMaps, groundTruth);
% fprintf('\nPrecision=%f Recall=%f MAE=%f\n', prec, rec, mae);  %Precision=%f Recall=%f 

OptData.frames = framesdata;
OptData.flow   = flowdata; %flowOrgdata,flowdata
OptData.nframe = numImages;
finalSpaTempOptSacyMaps = SpaTempSacyOptimization(frameSpaTempConSacyMaps,OptData,outPath);
% % Load for evaluation
% EstSpaTempSacyMaps = load(fullfile(outPath, 'SpaTempSacy-002-002-5-20-Alpa25', 'finalSpaTempOptSacyMaps.mat'));
% finalSpaTempOptSacyMaps = EstSpaTempSacyMaps.finalSpaTempOptSacyMaps;
    
% Evaluate SacyMaps
[prec, rec, mae] = evaluate_PR_MAE(finalSpaTempOptSacyMaps, groundTruth);
fprintf('\nfinal-SpaTempOptSacyMap: Precision=%f Recall=%f MAE=%f\n', prec, rec, mae);  

% Save Results
SpaTempOptSacyMapsSav = sprintf('%s/finalSpaTempOptSacyMaps.mat',outPath);
save(SpaTempOptSacyMapsSav,'finalSpaTempOptSacyMaps');


% %% Spatio-temperal Segmentation
% mode_MY = 'videoRapidSegment_Y';
% STSacyMapsLd = 'SP150FuseIou075-OpNo2-Ra05-SpTpOpt_Shen-Fnl';
% TData.imgs = framesdata;
% TData.flow = flowdata;
% TData.Orgflow = flowOrgdata;
% TData.flowmag = magnitude;
% % % Load framePixelList if it exist
% % framePixelListSav = fullfile(outPath, STSacyMapsLd, 'SpixelList.mat');
% % if exist(framePixelListSav, 'file')
% %     SPListLabels = load(framePixelListSav);
% %     framePixelList = SPListLabels.framePixelList;
% % end
% TData.superpixels = framePixelList; 
% 
% % % Load frameSpaTempConSacyMaps if it exist
% % SpaTempSacyMapsSav = fullfile(outPath, STSacyMapsLd, 'frameSpaTempSacyMaps.mat');
% % if exist(SpaTempSacyMapsSav, 'file')
% %     EstSpaTempSacyMaps = load(SpaTempSacyMapsSav);
% %     frameSpaTempConSacyMaps = EstSpaTempSacyMaps.frameSpaTempConSacyMaps;
% % end
% TData.saliency = frameSpaTempConSacyMaps;   %[frameSpaTempConSacyMaps, finalSpaTempOptSacyMaps;]
% % Compute Spatio-temperal object segments per frame
% SpaTempSegmentations = saliency2segmentation(TData, mode_MY);
% % Save the object segments per frame
% for idx = 1: length(SpaTempSegmentations)
%     SpaTempSeg = SpaTempSegmentations{idx};
%     outObjectPath = sprintf('%s/SpaTempObjSegsFram%03d.jpg',outPath,idx);    
%     imwrite(SpaTempSeg, outObjectPath);
% end
% 
% % % Load frameSpaTempConSacyMaps if it exist
% % finalSpaTempOptSacyMapsSav = fullfile(outPath, STSacyMapsLd, 'finalSpaTempOptSacyMaps.mat');
% % if exist(finalSpaTempOptSacyMapsSav, 'file')
% %     EstSpaTempSacyMaps = load(finalSpaTempOptSacyMapsSav);
% %     finalSpaTempOptSacyMaps = EstSpaTempSacyMaps.finalSpaTempOptSacyMaps;
% % end
% clear TData.saliency
% TData.saliency = finalSpaTempOptSacyMaps;
% % Compute the optimized Spatio-temperal object segments per frame
% SpaTempOptSegmentations = saliency2segmentation(TData, mode_MY);
% % Save the object segments per frame
% for idx = 1: length(SpaTempOptSegmentations)
%     SpaTempOptSeg = SpaTempOptSegmentations{idx};
%     outObjectPath = sprintf('%s/SpaTempOptObjSegsFram%03d.jpg',outPath,idx);    
%     imwrite(SpaTempOptSeg, outObjectPath);
% end


% %% Bounding Box Evaluation
% destGT = fullfile(imagePath, SeqName, 'gt');
% ground_bbox = gndtruth_extract(destGT);
% 
% % 1)SpaTempObjSegs --> Before SpaTempOptimization 
% detc_bbox = DetectBoxs(SpaTempSegmentations, ground_bbox, h, w, outPath, 15, 'SpaTempSegsBoxLabFram'); % BoundingBox extended length
% [ACC_one2all, ACC_one2one] = acc_eva(detc_bbox, ground_bbox);
% fprintf('\nExtend15ACC_one2all=%3.3f\n', ACC_one2all); 
% fprintf('\nExtend15ACC_one2all=%3.3f\n', ACC_one2one); 
% % 2)SpaTempOptObjSegs --> After SpaTempOptimization
% detc_bbox = DetectBoxs(SpaTempOptSegmentations, ground_bbox, h, w, outPath, 15, 'SpaTempOptSegsBoxLabFram'); % BoundingBox extended length
% [ACC_one2all, ACC_one2one] = acc_eva(detc_bbox, ground_bbox);
% fprintf('\nExtend15ACC_one2all=%3.3f\n', ACC_one2all); 
% fprintf('\nExtend15ACC_one2all=%3.3f\n', ACC_one2one); 
% 
% clear all
