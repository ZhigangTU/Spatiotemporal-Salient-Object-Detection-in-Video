function detc_bbox = DetectBoxs(SpaTempOptSegmentations, ground_bbox, h, w, outPath, EdL, SegSaveName)

detc_bbox = cell(size(SpaTempOptSegmentations));
for i = 1 : length(SpaTempOptSegmentations)
    stats = regionprops(imdilate(SpaTempOptSegmentations{i},strel('diamond',5)),'BoundingBox');  %strel('diamond',1);strel('disk',3,2)
    candidates = cell2mat(struct2cell(stats)');    
    if isempty(candidates) == 1
        candidates = [w/3, h/3, w/5, h/5];
    end
    indi = find(candidates(:,3)<10 | candidates(:,4)<10);
    candidates(indi,:) = [];
    % Expand the bounding box of the candidates
    Tempcandidates = [];
    for idx = 1 : size(candidates,1)
        Temp = candidates(idx,:);
        Tempcandidates = [Tempcandidates; [max(Temp(1)-EdL,1), max(Temp(2)-EdL,1), min(Temp(3)+(2*EdL),w), min(Temp(4)+(2*EdL),h)]];
    end
    candidates = Tempcandidates;       % struct('BoundingBox',Tempcandidates);
    
    if ~isempty(candidates)
        with_BB = insertShape(single(SpaTempOptSegmentations{i}),'Rectangle',candidates,'Color','green');
        with_BB = insertShape(with_BB,'Rectangle',ground_bbox{i},'Color','red');
        outObjectPath = sprintf('%s/%s%03d.jpg',outPath,SegSaveName,i);  
        imwrite(with_BB,outObjectPath)
    else
        outObjectPath = sprintf('%s/%s%03d.jpg',outPath,SegSaveName,i);  
        imwrite(SpaTempOptSegmentations{i},outObjectPath)
    end
    detc_bbox{i} = candidates;
end
