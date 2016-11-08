% function segmentation = Test_FastObjSeg(flowdata,magnitude,framesdata,h,w,numImages,outPath)

% addpath( genpath( '.' ) )

% foldername = fileparts( mfilename( 'fullpath' ) );

% The folder where the frames are stored in. Frames should be .jpg files
% and their names should be 8-digit sequential numbers (e.g. 00000001.jpg, 00000002.jpg etc)

currentPath = cd;

% input path
userName  = '001';
imagePath = fullfile('E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images','Riding-Horse');  % JHMDB\images\climb_stairs; data\input\Diving-Side
flowPath  = fullfile('E:\Action Recognition\Pose-CNN\cache\OFUCFSport','Riding-Horse');         % flowdata\Diving-Side
imageDir  = fullfile(imagePath, userName);
flowDir   = fullfile(flowPath, userName);
outPath   = fullfile(currentPath,'data\output\Riding-Horse',userName); % JHMDB\brush_hair

options.infolder = imageDir;

% The folder where all the outputs will be stored.
options.outfolder = outPath;

% The optical flow method to be used. Valid names are:
%   broxPAMI2011:     CPU-based optical flow.
%   sundaramECCV2010: GPU-based optical flow. Requires CUDA 5.0
%   options.flowmethod = 'broxPAMI2011';

% The superpixel oversegmentation method to be used. Valid names are:
%   Turbopixels
%   SLIC
options.superpixels = 'SLIC';

options.visualise = true;

% Print status messages on screen
options.vocal = true;

% options.ranges:
%   A matlab array of length S+1, containing the number for the 
%   first frame of each shot (where S is the total count of shots
%   inside the options.infolder). 
%   The last element of the array should be equal to the total 
%   number of frames + 1.
%   options.ranges = [ 1, 33, 48 ];

% options.positiveRanges:
% A matlab array containing the shots to be processed
options.positiveRanges = [ 1, 2 ];

% If the frames are larger than options.maxedge in either height or width, they will be resized to fit a (maxedge x maxedge) window. 
% This greatly decreases the optical flow computation cost, without (typically) degrading the segmentation accuracy too much. 
% If resizing the frame is not desirable, set options.maxedge = inf
options.maxedge = 1024;

% Use default params. For specific value info check inside the function
params = getDefaultParams();

% Create folder to save the segmentation
segmfolder = fullfile(options.outfolder, 'segmentations', 'VideoRapidSegment' );
if( ~exist( segmfolder, 'dir' ) )
    mkdir( segmfolder )
end;

shot = 1;

%% Get training images
[fileNames, numImages] = gzget_training_images(imagePath, userName);

%% read first image to get the size info
testImage = imread(fileNames{1});
if size(testImage, 3) > 1
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

%% read image file, to establish the data matrix
framesdata = cell(numImages, 1);
% framesName = cell(numImages, 1);
for fileIndex = 1:numImages
    currentImage = imread(fileNames{fileIndex});
    framesdata{fileIndex} = currentImage;
end

%% read flow data, u v and norm %%
flowOrgdata  = cell(1, numImages-1);
flowdata  = cell(1, numImages-1);
magnitude = cell(1, numImages-1);
 
for j = 1:numImages - 1
    infName0 = sprintf('flowu%d.mat',j);
    infName00 = fullfile(flowDir, infName0);
    load(infName00,'u');
    uflow = int16(u);
    
    infName1 = sprintf('flowv%d.mat',j);
    infName11 = fullfile(flowDir, infName1);
    load(infName11,'v');
    vflow = int16(v);
    
    flowOrgdata{j} = cat(3,u,v);
    flowdata{j}  = cat(3,uflow,vflow);
    Flowgradient = getFlowGradient(flowdata{j});
    magnitude{j} = getMagnitude(Flowgradient);
end

data.Orgflow = flowOrgdata;
data.flow = flowdata;
data.imgs = framesdata;
data.flowmag = magnitude;

% Load superpixels (or compute if not found)
data.superpixels = loadSuperpixels(options, shot);
if(isempty(data.superpixels))
    data.superpixels = computeSuperpixels(options, shot, numImages, framesdata);
end

data.id = shot;
segmentationCom = videoRapidSegmentTu(options, params, data);
segmentationOrg = segmentationCom(:,:,1);
segmentationUin = segmentationCom(:,:,2);
segmentation = cell(length(segmentationOrg),1);
for idx = 1: length(segmentationOrg)    
    segOrg = segmentationOrg{idx};
    segUin = segmentationUin{idx};
    [m1, n1] = size(segOrg);
    segRef = segOrg;
    [LOrg,numtargetOrg]=bwlabel(segOrg,8);
    [LUin,numtargetUin]=bwlabel(segUin,8);
    if numtargetOrg < 1   % No objects are detected in OrgFlow
        segRef = segUin;
    else                  % Refine OrgFlow objects with UinFlow objects
        for iOTarg = 1:numtargetOrg
            [rtemp,ctemp] = find(LOrg==iOTarg);
            minH=min(rtemp);
            minW=min(ctemp);
            maxH=max(rtemp);
            maxW=max(ctemp);
            minH = max(1,minH);
            minW = max(1,minW);
            maxH = min(m1,maxH);
            maxW = min(n1,maxW);
            % TargetOBox = [minH minW maxH maxW];
            box1 = [minW,minH,(maxW-minW),(maxH-minH)];   %[X,Y,WIDTH,HEIGHT]
            for iUTarg = 1:numtargetUin
                [Urtemp,Uctemp] = find(LUin==iUTarg);
                UminH=min(Urtemp);
                UminW=min(Uctemp);
                UmaxH=max(Urtemp);
                UmaxW=max(Uctemp);
                UminH = max(1,UminH);
                UminW = max(1,UminW);
                UmaxH = min(m1,UmaxH);
                UmaxW = min(n1,UmaxW);
                box2 = [UminW,UminH,(UmaxW-UminW),(UmaxH-UminH)];
                iou = inters_union(box1,box2);
                if iou > 0
                    RminH = min(minH,UminH);
                    RminW = min(minW,UminW);
                    RmaxH = max(maxH,UmaxH);
                    RmaxW = max(maxW,UmaxW);
                    segRef(RminH:RmaxH,RminW:RmaxW) = segOrg(RminH:RmaxH,RminW:RmaxW)+segUin(RminH:RmaxH,RminW:RmaxW);
                end
            end            
        end
    end
    
    segTargets = segRef>0;  
    
    % Delect and clean noisy pixels
    [L,numtarget]=bwlabel(segTargets,8);  
    stats = regionprops(L,'Area');     % Compute the size of each connected region
    area  = cat(1,stats.Area);     
    index = find(area > 8^2);          % Finding the index of the region that its size larger than a threshold (e.g.10^2)   
    % Obtaining the indexed regions which eliminate small noisy regions 
    segTargetsC = ismember(L,index(:));    
    segmentation{idx} = segTargetsC;    
end

for idx = 1: length(segmentation)
    seg = segmentation{idx};
    outObjectPath = sprintf('%s/FastObjSeg-Fram%03d.jpg',outPath, idx);    
    imwrite(seg, outObjectPath);
end

% Save output
filename = fullfile( segmfolder, sprintf('segmentation.mat'));
save(filename, 'segmentation', '-v7.3');   

% rmpath( genpath( '.' ) )

clear all
