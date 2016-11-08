param=[];
param.lhandposition=13; % pose joints positions in the structure (JHMDB pose format)
param.rhandposition=12;
param.upbodypositions=[1 2 3 4 5 6 7 8 9 12 13];
param.lside = 40; % length of part box side (also depends on the human scale)

param.impath = 'E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images\Walk-Front'; % input images path (one folder per video); 'Data/JHMDB/images'; 'Data/UCF-Sports/images';clean/final;
param.imext  = '.png';                % input image extension type

% get video names
video_names = dir(param.impath);
video_names = {video_names.name};
video_names = video_names(~ismember(video_names,{'.','..'}));

% compute_OF(video_names,param);  % compute optical flow between adjacent frames

% dname=sprintf('%s/%s','C:/AR Dataset/Data/UCF101','OF');
% if ~exist(dname,'dir')
%     mkdir(dname);
% end

fprintf('\n------ Object Detection ------\n')

nb_vid = length(video_names);
impath = param.impath;
OFpath = param.ofpath;
Sapath = param.cachepath;

for vi = 1:nb_vid       % parfor/for
    vidname = video_names{vi};
    vidnameTep = '001';
    if strcmp(vidname,vidnameTep)
        continue;
    end
    fprintf('\nObject Detection on the video sequence: %s\n',vidname)
    
    Test_TRPCASacyMAP(vidname)
    
end