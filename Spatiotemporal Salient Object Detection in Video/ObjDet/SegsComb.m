function segmentation = SegsComb(segmentationCom)

if size(segmentationCom,3)>1
    % Refine the two-stream segments
    segmentationUin = segmentationCom(:,:,1);
    segmentationOrg = segmentationCom(:,:,2);
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
else
    segmentation = segmentationCom;
end
