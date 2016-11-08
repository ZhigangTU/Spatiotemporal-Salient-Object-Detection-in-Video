function TargetObjects = TargetRPCAMotionFuse(frameRPCABoxs, frameMotionBoxs, nframes, frames, magflows,outPath)

frameRPCABoxsDel = cell(nframes, 1);
frameRPCAReBoxs = cell(nframes, 1);

Height = size(frames{1},1);
Width  = size(frames{1},2);

for j = 1:nframes
    fRPCABoxsDel = [];
    fRPCABoxsDelR = [];
    
    fRPCABox = frameRPCABoxs{j};
    fMotionBox = frameMotionBoxs{j};    
    iRPCARow = size(fRPCABox,1);
    iMotionRow = size(fMotionBox,1);
    
%     iframeRPCA = frames{j};
%     for i = 1:iRPCARow
%             coortemp = fRPCABox(i,:);
%             hmin = coortemp(2);
%             wmin = coortemp(3);
%             hmax = coortemp(4);
%             wmax = coortemp(5);
%             % Labeling the boxes
%             iframeRPCA(hmin : hmin, wmin : wmax) = 255;
%             iframeRPCA(hmax : hmax, wmin : wmax) = 255;
%             iframeRPCA(hmin : hmax, wmin : wmin) = 255;
%             iframeRPCA(hmin : hmax, wmax : wmax) = 255;
%      end
%      figure; imshow(uint8(iframeRPCA))
%      iframeMotion = frames{j};
%      for i = 1:iMotionRow
%             coortemp = fMotionBox(i,:);
%             hmin = coortemp(2);
%             wmin = coortemp(3);
%             hmax = coortemp(4);
%             wmax = coortemp(5);
%             % Labeling the boxes
%             iframeMotion(hmin : hmin, wmin : wmax) = 255;
%             iframeMotion(hmax : hmax, wmin : wmax) = 255;
%             iframeMotion(hmin : hmax, wmin : wmin) = 255;
%             iframeMotion(hmin : hmax, wmax : wmax) = 255;
%      end
%      figure; imshow(uint8(iframeMotion))
     

    AreaBoxs = [];
    MaxRPCABox = [];
    for iiR = 1:iRPCARow
        iiRPCABox = fRPCABox(iiR,:);
        Areatemp = (iiRPCABox(5)-iiRPCABox(3))*(iiRPCABox(4)-iiRPCABox(2));
        AreaBoxs  = [AreaBoxs;Areatemp];
    end
    iarMax = find(AreaBoxs==max(AreaBoxs));     % The index of the smallest box
    iarMax = iarMax(1);
    MaxRPCABox = fRPCABox(iarMax,:);
    if j < nframes
        magflow = magflows{j};
    else 
        magflow = magflows{nframes-1};
    end
    Aremagflow = magflow(MaxRPCABox(2):MaxRPCABox(4),MaxRPCABox(3):MaxRPCABox(5));      
    MxmagFlowAre = max(Aremagflow(:));
    MemagFlowAre = mean(Aremagflow(:));
        
    % Delete the RPCAboxes that have no iou with MotionBoxes
    if (MxmagFlowAre<=2.0) && (MemagFlowAre<=0.020) && (AreaBoxs(iarMax)>=(1/6)*Height*Width)   % Human with small walking speed
        fRPCABoxsDelR = [fRPCABoxsDelR; MaxRPCABox(2) MaxRPCABox(3) MaxRPCABox(4) MaxRPCABox(5)];
    else
        for iR = 1:iRPCARow
            ifRPCABox = fRPCABox(iR,:);
            box1 = [ifRPCABox(3),ifRPCABox(2),(ifRPCABox(5)-ifRPCABox(3)),(ifRPCABox(4)-ifRPCABox(2))];   %[X,Y,WIDTH,HEIGHT]
            nInside = 0;
            for iM = 1:iMotionRow
                ifMotionBox = fMotionBox(iM,:);           
                box2 = [ifMotionBox(3),ifMotionBox(2),(ifMotionBox(5)-ifMotionBox(3)),(ifMotionBox(4)-ifMotionBox(2))];
                iou = inters_union(box1,box2);
                if iou > 0 && nInside < 1      % >0.1
                    fRPCABoxsDelR = [fRPCABoxsDelR;ifRPCABox(2) ifRPCABox(3) ifRPCABox(4) ifRPCABox(5)];   
                    nInside = nInside+1;
                end
            end
        end
    end
    
    if isempty(fRPCABoxsDelR)
%         AreaBoxs = [];
%         for iiR = 1:iRPCARow
%             iiRPCABox = fRPCABox(iiR,:);
%             Areatemp = (iiRPCABox(5)-iiRPCABox(3))*(iiRPCABox(4)-iiRPCABox(2));
%             AreaBoxs  = [AreaBoxs;Areatemp];
%         end
%         iarMax = find(AreaBoxs==max(AreaBoxs));     % The index of the smallest box 
%         MaxRPCABox = fRPCABox(iarMax,:);
%         if j < nframes
%             magflow = magflows{j};
%         else 
%             magflow = magflows{nframes-1};
%         end
%         Aremagflow = magflow(MaxRPCABox(2):MaxRPCABox(4),MaxRPCABox(3):MaxRPCABox(5));      
%         MxmagFlowAre = max(Aremagflow(:));
%         MemagFlowAre = mean(Aremagflow(:));
        if (MxmagFlowAre<=2.0) && (MemagFlowAre<=0.020)    % Human with small walking speed
            fRPCABoxsDelR = [fRPCABoxsDelR; MaxRPCABox(2) MaxRPCABox(3) MaxRPCABox(4) MaxRPCABox(5)];
        end
        for i = 1:iMotionRow
            coortemp = fMotionBox(i,:);
            fRPCABoxsDelR = [fRPCABoxsDelR; coortemp(2) coortemp(3) coortemp(4) coortemp(5)];
        end
    end
    
    fRPCABoxsDelR = gzDeleteInsideblock(fRPCABoxsDelR);
    
    if isempty(fRPCABoxsDelR)
        for i = 1:iMotionRow
            coortemp = fMotionBox(i,:);
            fRPCABoxsDel = [fRPCABoxsDel; coortemp(1) coortemp(2) coortemp(3) coortemp(4) coortemp(5)];
        end
    else
        iRPCADelRow = size(fRPCABoxsDelR,1);
        for i = 1:iRPCADelRow
            coortemp = fRPCABoxsDelR(i,:);
            fRPCABoxsDel = [fRPCABoxsDel; j coortemp];
        end       
    end
    
    % Boxes combination refinement 
    TempImg = zeros(Height, Width);
    for ii = 1:size(fRPCABoxsDel,1)
        boxtemp = fRPCABoxsDel(ii,:);
        TempImg(boxtemp(2):boxtemp(4),boxtemp(3):boxtemp(5))=1;
    end
    EdgImg = (TempImg>0);
    [L,numtarget]=bwlabel(EdgImg,8);
    fRPCABoxsDel = [];
    for nT = 1:numtarget
        [rtemp,ctemp] = find(L==nT);
        minH=min(rtemp);
        minW=min(ctemp);
        maxH=max(rtemp);
        maxW=max(ctemp);

        minH = max(1,minH);
        minW = max(1,minW);
        maxH = min(Height,maxH);
        maxW = min(Width,maxW);

        fRPCABoxsDel = [fRPCABoxsDel; j minH minW maxH maxW];
    end
    
    iframeRPCADel = frames{j};
    iRPCADelRow = size(fRPCABoxsDel,1);
    for i = 1:iRPCADelRow
        coortemp = fRPCABoxsDel(i,:);
        hmin = coortemp(2);
        wmin = coortemp(3);
        hmax = coortemp(4);
        wmax = coortemp(5);
        % Labeling the boxes
        iframeRPCADel(hmin : hmin, wmin : wmax) = 255;
        iframeRPCADel(hmax : hmax, wmin : wmax) = 255;
        iframeRPCADel(hmin : hmax, wmin : wmin) = 255;
        iframeRPCADel(hmin : hmax, wmax : wmax) = 255;
    end
%     figure; imshow(uint8(iframeRPCADel))
    outObjectPath = sprintf('%s/frameRPCABox%05d.jpg',outPath, j);
    imwrite(uint8(iframeRPCADel), outObjectPath);
           
    frameRPCABoxsDel{j} = fRPCABoxsDel;    
end

for j = 1:nframes
    fRPCABoxs = [];
    fRPCAMBoxs = [];
    
    fRPCABoxRef = frameRPCABoxsDel{j};
    fMotionBox = frameMotionBoxs{j};    
    iRPCARow = size(fRPCABoxRef,1);
    iMotionRow = size(fMotionBox,1);
    
    for iR = 1:iRPCARow
        iRPCABox = fRPCABoxRef(iR,:);
        box1 = [iRPCABox(3),iRPCABox(2),(iRPCABox(5)-iRPCABox(3)),(iRPCABox(4)-iRPCABox(2))];   %[X,Y,WIDTH,HEIGHT]
        ar1 = box1(:,3).*box1(:,4);
        nIoU = 0;
        nInside = 0;
        for iM = 1:iMotionRow
            coortempMo = fMotionBox(iM,1:5);
            box2 = [coortempMo(3),coortempMo(2),(coortempMo(5)-coortempMo(3)),(coortempMo(4)-coortempMo(2))];
            ar2 = box2(:,3).*box2(:,4);
            iou = inters_union(box1,box2);
            if (iou>0 && iou<1)
                nIoU = nIoU+1;
                if ar1<=ar2
                    if ar2<(Height*Width*0.5)
                        if size(fRPCABoxs,1)>=1
                            boxtempM = repmat(coortempMo(:)', [size(fRPCABoxs,1) 1]);
                            boxtempD = abs(boxtempM - fRPCABoxs);
                            boxtempD = sum(boxtempD,2);
                            boxtempDLab = (boxtempD>0);
                            if sum(boxtempDLab(:))==(size(fRPCABoxs,1))
                                fRPCABoxs = [fRPCABoxs;coortempMo];
                            end
                        else
                            fRPCABoxs = [fRPCABoxs;coortempMo];
                        end
                    else
                        if size(fRPCABoxs,1)>=1
                            boxtempM = repmat(coortempMo(:)', [size(fRPCABoxs,1) 1]);
                            boxtempD = abs(boxtempM - fRPCABoxs);
                            boxtempD = sum(boxtempD,2);
                            boxtempDLab = (boxtempD>0);
                            if sum(boxtempDLab(:))==(size(fRPCABoxs,1))
                                hmin = min(coortempMo(2)+round(1/6*(coortempMo(4)-coortempMo(2))),Height-2);
                                wmin = min(coortempMo(3)+round(1/6*(coortempMo(5)-coortempMo(3))),Width-2);
                                hmax = max(coortempMo(4)-round(1/6*(coortempMo(4)-coortempMo(2))),hmin+2);
                                wmax = max(coortempMo(5)-round(1/6*(coortempMo(5)-coortempMo(3))),wmin+2);
                                coordMo = [hmin wmin hmax wmax];
                                fRPCABoxs = [fRPCABoxs;j coordMo];
                            end
                        else
                            hmin = min(coortempMo(2)+round(1/6*(coortempMo(4)-coortempMo(2))),Height-2);
                            wmin = min(coortempMo(3)+round(1/6*(coortempMo(5)-coortempMo(3))),Width-2);
                            hmax = max(coortempMo(4)-round(1/6*(coortempMo(4)-coortempMo(2))),hmin+2);
                            wmax = max(coortempMo(5)-round(1/6*(coortempMo(5)-coortempMo(3))),wmin+2);
                            coordMo = [hmin wmin hmax wmax];
                            fRPCABoxs = [fRPCABoxs;j coordMo];
                        end
                    end    
                elseif (ar1>=3*ar2 && ar1>=(Height*Width*0.4))   %RPCA too large, larger the small Motion box
                    hmin = max(coortempMo(2)-round(1/2*(coortempMo(4)-coortempMo(2))),1);
                    wmin = max(coortempMo(3)-round(1/2*(coortempMo(5)-coortempMo(3))),1);
                    hmax = min(coortempMo(4)+round(1/2*(coortempMo(4)-coortempMo(2))),Height);
                    wmax = min(coortempMo(5)+round(1/2*(coortempMo(5)-coortempMo(3))),Width);
                    coordMo = [hmin wmin hmax wmax];
                    fRPCABoxs = [fRPCABoxs;j coordMo];
                else
                    minH = min(iRPCABox(2),coortempMo(2));
                    minW = min(iRPCABox(3),coortempMo(3));
                    maxH = max(iRPCABox(4),coortempMo(4));
                    maxW = max(iRPCABox(5),coortempMo(5));
                    boxcom = [minH minW maxH maxW];
                    fRPCABoxs = [fRPCABoxs;j boxcom];
                end
            elseif nIoU>0
                continue;            
            elseif nIoU==0 && nInside<1
                fRPCABoxs = [fRPCABoxs;iRPCABox(1) iRPCABox(2) iRPCABox(3) iRPCABox(4) iRPCABox(5)];
                nInside = nInside+1;
            end
        end  
        % Delete the inside boxes
        fRPCABoxsP = fRPCABoxs(:,2:5);
        % Check the same box       
        if size(fRPCABoxsP,1)>=2
            fRPCABoxsPD = [];
            iBoxsPRow = size(fRPCABoxsP,1);
            iLBoxs = fRPCABoxsP(iBoxsPRow,:);
            for iBR = 1:iBoxsPRow-1
                BoxD = fRPCABoxsP(iBR,:)-iLBoxs;
                if sum(BoxD(:))==0
                    fRPCABoxsP(iBR,:) = [0 0 0 0];
                end
            end
            for iBRD = 1:iBoxsPRow
                iBoxRD = fRPCABoxsP(iBRD,:);
                if sum(iBoxRD(:))>0
                    fRPCABoxsPD = [fRPCABoxsPD;iBoxRD];
                end
            end
            fRPCABoxsP = fRPCABoxsPD;
        end        
        fRPCABoxsP = gzDeleteInsideblock(fRPCABoxsP);
        if isempty(fRPCABoxsP)
            fRPCABoxsP = fRPCABoxs(1,2:5);
        end
        fRPCABoxs = [j*ones(size(fRPCABoxsP,1),1) fRPCABoxsP];
        % Delete the small box belongs to RPCA boxes
        BoxsRef = [];
        AreaRef = [];
        nSmallRPCAB = 0;
        for ii = 1:size(fRPCABoxs)
            iiBox = fRPCABoxs(ii,:);
            Boxstemp = [iiBox(3),iiBox(2),(iiBox(5)-iiBox(3)),(iiBox(4)-iiBox(2))];
            BoxsRef  = [BoxsRef;Boxstemp];
            Areatemp = Boxstemp(:,3).*Boxstemp(:,4);
            AreaRef  = [AreaRef;Areatemp];
        end
        iarMin = find(AreaRef==min(AreaRef));     % The index of the smallest box
        MinCoord = fRPCABoxs(iarMin,:);           % The coordinates of the smallest box
        MinBox = BoxsRef(iarMin,:);
        MinCoord = MinCoord(1,:);
        MinBox = MinBox(1,:);
        for ii = 1:size(fRPCABoxs)
            iiBox = BoxsRef(ii,:);
            iou = inters_union(MinBox,iiBox);
            if (iou>0 && iou<1)
                nSmallRPCAB = nSmallRPCAB+1;
            end
        end       
        if (sum(MinCoord-iRPCABox)==0) && (nSmallRPCAB>0)  % Whether MinCoord belongs to RPCA boxes
            fRPCABoxs(iarMin,:) = [];      
        end
    end
    fRPCABoxsP2 = fRPCABoxs(:,2:5);
    fRPCABoxsP2 = gzDeleteInsideblock(fRPCABoxsP2);
    if isempty(fRPCABoxsP2)
        fRPCABoxsP2 = fRPCABoxs(1,2:5);
    end
    fRPCABoxs = [j*ones(size(fRPCABoxsP2,1),1) fRPCABoxsP2];
    
    for ii = 1:iMotionRow
        boxtemp = fMotionBox(ii,1:5);
        box1 = [boxtemp(3),boxtemp(2),(boxtemp(5)-boxtemp(3)),(boxtemp(4)-boxtemp(2))];   %[X,Y,WIDTH,HEIGHT]
        for jj = 1:size(fRPCABoxs,1)               
            boxcmp = fRPCABoxs(jj,:);                    
            box2 = [boxcmp(3),boxcmp(2),(boxcmp(5)-boxcmp(3)),(boxcmp(4)-boxcmp(2))];
            iou = inters_union(box1,box2);
            if iou > 0
                minH   = min(boxtemp(2),boxcmp(2));
                minW   = min(boxtemp(3),boxcmp(3));
                maxH   = max(boxtemp(4),boxcmp(4));
                maxW   = max(boxtemp(5),boxcmp(5));
                boxcom = [minH minW maxH maxW];
                fRPCAMBoxs = [fRPCAMBoxs;j boxcom];
                boxtemp = [j boxcom];
                box1    = [boxcom(2),boxcom(1),(boxcom(4)-boxcom(2)),(boxcom(3)-boxcom(1))];
            end
        end
        if isempty(fRPCAMBoxs)
            continue;
        else
            fRPCAMBoxsP4 = fRPCAMBoxs(:,2:5);
            fRPCAMBoxsP4 = gzDeleteInsideblock(fRPCAMBoxsP4);
            if isempty(fRPCAMBoxsP4)
                fRPCAMBoxsP4 = fRPCAMBoxs(1,2:5);
            end
            fRPCAMBoxs   = [j*ones(size(fRPCAMBoxsP4,1),1) fRPCAMBoxsP4];
            
            boxtempM = repmat(fRPCAMBoxs(:)', [size(fRPCABoxs,1) 1]);
            boxtempD = abs(boxtempM - fRPCABoxs);
            boxtempD = sum(boxtempD,2);
            boxtempDLab = (boxtempD>0);
            if sum(boxtempDLab(:))==(size(fRPCABoxs,1))
                fRPCABoxs = [fRPCABoxs;fRPCAMBoxs];
            end
            
            fRPCABoxsP3  = fRPCABoxs(:,2:5);
            fRPCABoxsP3  = gzDeleteInsideblock(fRPCABoxsP3);
            if isempty(fRPCABoxsP3)
                fRPCABoxsP3 = fRPCABoxs(1,2:5);
            end
            fRPCABoxs    = [j*ones(size(fRPCABoxsP3,1),1) fRPCABoxsP3];
            fRPCAMBoxs   = [];
            fRPCABoxsP3  = [];
        end
    end

    % Boxes combination refinement 
    TempImg = zeros(Height, Width);
    for ii = 1:size(fRPCABoxs,1)
        boxtemp = fRPCABoxs(ii,:);
        TempImg(boxtemp(2):boxtemp(4),boxtemp(3):boxtemp(5))=1;
    end
    EdgImg = (TempImg>0);
    [L,numtarget]=bwlabel(EdgImg,8);
    fRPCABoxs = [];
    for nT = 1:numtarget
        [rtemp,ctemp] = find(L==nT);
        minH=min(rtemp);
        minW=min(ctemp);
        maxH=max(rtemp);
        maxW=max(ctemp);

        minH = max(1,minH);
        minW = max(1,minW);
        maxH = min(Height,maxH);
        maxW = min(Width,maxW);
        if (maxH-minH)*(maxW-minW) >= 0.5*Height*Width
            BoxH = maxH-minH;
            BoxW = maxW-minW;            
            minH = min(minH+round(1/6*BoxH),Height-2);
            minW = min(minW+round(1/6*BoxW),Width-2);
            maxH = max(maxH-round(1/6*BoxH),minH+2);
            maxW = max(maxW-round(1/6*BoxW),minW+2);
            fRPCABoxs = [fRPCABoxs; j minH minW maxH maxW];
        else
            fRPCABoxs = [fRPCABoxs; j minH minW maxH maxW];
        end
    end
        
    
    iframeRPCARef = frames{j};
    iRPCARefRow = size(fRPCABoxs,1);
    for i = 1:iRPCARefRow
        coortemp = fRPCABoxs(i,:);
        hmin = coortemp(2);
        wmin = coortemp(3);
        hmax = coortemp(4);
        wmax = coortemp(5);
        % Labeling the boxes
        iframeRPCARef(hmin : hmin, wmin : wmax) = 255;
        iframeRPCARef(hmax : hmax, wmin : wmax) = 255;
        iframeRPCARef(hmin : hmax, wmin : wmin) = 255;
        iframeRPCARef(hmin : hmax, wmax : wmax) = 255;
    end
    figure; imshow(uint8(iframeRPCARef))
    outObjectPath = sprintf('%s/frameRPCAMotionRef%05d.jpg',outPath, j);
    imwrite(uint8(iframeRPCARef), outObjectPath);
    
    frameRPCAReBoxs{j} = fRPCABoxs;
end

TargetObjects = frameRPCAReBoxs;

% ObjectsCoord = sprintf('%s/frameRPCADel.mat',outPath);
% save(ObjectsCoord,'TargetObjects');

