function segmentation = saliency2segmentation(data,mode_MY)

%method of averaging saliency 'mean' or 'median'
ave_mod = 'mean';
%select Segment method : videoRapidSegment_Y/videoRapidSegment_M
seg_method = mode_MY;

options.superpixels = 'SLIC';
options.visualise = false;
options.vocal = true;
options.maxedge = max(size(data.imgs{1}));
params = getDefaultParams();

if strcmp(mode_MY,'videoRapidSegment_Y')
    params.appearanceUnaryWeight = 0.35; %0.35
    params.locationWeight = 0.5; %0.8
    params.temporalWeight = 50; %200
    params.spatialWeight = 500; %400
else
    params.appearanceUnaryWeight = 0.3; %0.1 very aggressive eliminate unwanted segmentations
                                        %0.3 usual setting
                                        %0.5 seems better
    params.locationWeight = 0.8; %2
                               %3 more balanced
    params.temporalWeight = 50;
    params.spatialWeight = 500;
end

data.id = 1;

new_superpixel = cell(size(data.saliency));
new_saliency = cell(size(data.saliency));
for i = 1 : length(data.superpixels)
    temp_sup = data.superpixels{i};
    temp_sal = data.saliency{i};
    dummy_sup = zeros(size(data.saliency{i}));
    dummy_sal = zeros(size(data.superpixels{i}));
    for j = 1 : length(temp_sup)
        dummy_sup(temp_sup{j}) = j;
        if strcmp(ave_mod,'mean')
            dummy_sal(j) = mean(temp_sal(temp_sup{j}));
        else
            dummy_sal(j) = median(temp_sal(temp_sup{j}));
        end
    end
    new_superpixel{i} = dummy_sup;
    new_saliency{i} = dummy_sal;
end


if strcmp(seg_method,'videoRapidSegment_Y')
    data.saliency = new_saliency;
    data.superpixels = cellfun(@(x) uint16(x),new_superpixel,'UniformOutput',false);
    segmentation = videoRapidSegment_Y( options, params, data );
else
    data.saliency = new_saliency;
    data.superpixels = cellfun(@(x) uint16(x),new_superpixel,'UniformOutput',false);
    segmentation = videoRapidSegment_M( options, params, data );
end

for idx = 1: size(segmentation,1)
    % Delect and clean noisy pixels
    frameSeg = segmentation{idx};
    [L,numtarget]=bwlabel(frameSeg,8);
    stats = regionprops(L,'Area');  % Compute the size of each connected region
    area = cat(1,stats.Area);     
    index = find(area > 20^2);       % Finding the index of the region that its size larger than a threshold (e.g.10*10)   
    % Obtaining the indexed regions which eliminate small noisy regions 
    frameSegClean = ismember(L,index(:));   
    if sum(frameSegClean(:)) < 1 
        if idx > 1
            segmentation{idx} = segmentation{idx-1};
        else
            frameSegLab = zeros(size(frameSegClean));
            [h,w] = size(frameSegLab);
            frameSegLab(round(h/3):round(2*h/3),round(w/3):round(2*w/3))=1;
            frameSegLab = frameSegLab>0;
            segmentation{idx}=frameSegLab;
        end
    else
        segmentation{idx} = frameSegClean;
    end
end

end
