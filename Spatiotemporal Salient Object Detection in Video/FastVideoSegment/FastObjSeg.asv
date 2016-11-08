function FOSsegmentation = FastObjSeg(data,numImages,flowPath,imagePath,outPathFOS,SeqName,segmfolder)

% input path
imageDir  = fullfile(imagePath, SeqName);
flowDir   = fullfile(flowPath, SeqName);
outPath   = outPathFOS; % JHMDB\brush_hair

options.infolder = imageDir;
options.SPoutfolder = flowPath;

% The superpixel oversegmentation method to be used. Valid names are:
%   Turbopixels
%   SLIC
options.superpixels = 'SLIC'; %'Turbopixels','SLIC';
options.maxedge = inf;
options.vocal = true;
options.visualise = true;
options.SeqName = SeqName;
options.flowmethod = 'broxPAMI2011'; %broxPAMI2011, NL, MDP, Brox

% A matlab array containing the shots to be processed
options.positiveRanges = [1, 2];

params = getSegTrackParams(SeqName); %getDefaultParams();getSegTrackParams; getSegTrackV2Params;

% data.Orgflow = flowOrgdata;
% data.flow = flowdata;
% data.imgs = framesdata;
% data.flowmag = magnitude;

shot = data.shot;

% Load superpixels (or compute if not found)
data.superpixels = loadSuperpixels(options, flowPath, shot);
if(isempty(data.superpixels))
    data.superpixels = computeSuperpixels(options, shot, numImages, data.imgs);
end

% Fast object segmentation 
data.id = shot;
segmentationCom = videoRapidSegmentTu(options, params, data);

% Combining Segments
FOSsegmentation = SegsComb(segmentationCom);
for idx = 1: length(FOSsegmentation)
    seg = FOSsegmentation{idx};
    outObjectPath = sprintf('%s/FastObjSeg-Fram%03d.jpg',outPath, idx);    
    imwrite(seg, outObjectPath);
end

segmentationLarg = getLargestSegmentAndNeighbours( FOSsegmentation );
for idx = 1: length(segmentationLarg)
    SpaTempOptSegLarg = segmentationLarg{idx};
    outObjectPath = sprintf('%s/FOSSegsLargFram%03d.jpg',outPath, idx);    
    imwrite(SpaTempOptSegLarg, outObjectPath);
end

% Save output
filename = fullfile( segmfolder, sprintf('segmentation.mat'));
save(filename, 'FOSsegmentation', '-v7.3');   
filename = fullfile( segmfolder, sprintf('segmentationLarg.mat'));
save(filename, 'segmentationLarg', '-v7.3'); 

% % Saliency evaluation on SegTrak dataset.
% avgMislabelled = getAverageMislabelledPixels( options, segmentationSeg1);
% fprintf( 'Average number of mislabelled pixels for %s: %i\n', SeqName, avgMislabelled ); % videoid{ shot }    
% avgMislabelled = getAverageMislabelledPixels( options, FOSsegmentation);
% fprintf( 'Average number of mislabelled pixels for %s: %i\n', SeqName, avgMislabelled ); % videoid{ shot } 
% avgMislabelled = getAverageMislabelledPixels( options, segmentationLarg);
% fprintf( 'Average number of mislabelled pixels for %s: %i\n', SeqName, avgMislabelled ); % videoid{ shot }

