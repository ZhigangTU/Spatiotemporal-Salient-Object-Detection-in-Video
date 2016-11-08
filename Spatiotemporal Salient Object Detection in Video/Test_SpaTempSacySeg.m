
param.impath = 'F:\Action-Recog Experiments\Pose-CNN\Data\UCFSports\images';

% get video names
video_names = dir(param.impath);
video_names = {video_names.name};
video_names = video_names(~ismember(video_names,{'.','..'}));

nb_vid = length(video_names);

% UCFSports Database
for vi = 1:nb_vid
    % get vedio sequences list in the current video
    vidname = video_names{vi};
    
    subimpath = sprintf('%s/%s',param.impath,vidname);
    subvideo_names = dir(subimpath);
    subvideo_names = {subvideo_names.name};
    subvideo_names = subvideo_names(~ismember(subvideo_names,{'.','..'}));
    
    subnb_vid = length(subvideo_names);
    
    for svi = 1:subnb_vid      
        SeqName = subvideo_names{svi};

        Test_TRPCASacyMAP(vidname, SeqName, svi)
    end    
end

% % SegTrack Database   
% for vi = 1:nb_vid 
%     vidname = video_names{vi};
%     Test_STSacySegTrack(vidname)
% end
    

