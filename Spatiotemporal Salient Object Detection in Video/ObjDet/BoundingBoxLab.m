function ObjectBox = BoundingBoxLab(GEdgeMap, j, nframe, Imframe, frameObjectBox, outPath, NumT)
if nargin < 7
    NumT = 1;
end

[L,numtarget]=bwlabel(GEdgeMap,8);
[m1, n1] = size(GEdgeMap);
TargetCoord = [];
TargetCoordR = [];
TargetNoOverlap = [];

if( ~exist(outPath, 'dir'))
    mkdir(outPath)
end

if numtarget>=1
    TargetTempSave = [];
    for nT = 1:numtarget,
        [rtemp,ctemp] = find(L==nT);
        minH=min(rtemp);
        minW=min(ctemp);
        maxH=max(rtemp);
        maxW=max(ctemp);
        
        if NumT == 1
            nExtend = 4;  % JHMDB-->2; Sports-->4
        else
            nExtend = 2;  % JHMDB-->1; Sports-->2
        end        
        minH = max(1,minH - nExtend);
        minW = max(1,minW - nExtend);
        maxH = min(m1,maxH + nExtend);
        maxW = min(n1,maxW + nExtend);
        TargetTempSave = [TargetTempSave; minH minW maxH maxW];
    end
    
    Row = size(TargetTempSave,1);
    Imframe0 = Imframe;
    for i = 1:Row
        coortemp = TargetTempSave(i,:);
        hmin = coortemp(1);
        wmin = coortemp(2);
        hmax = coortemp(3);
        wmax = coortemp(4);
        % Labeling the boxes
        Imframe0(hmin : hmin, wmin : wmax) = 255;
        Imframe0(hmax : hmax, wmin : wmax) = 255;
        Imframe0(hmin : hmax, wmin : wmin) = 255;
        Imframe0(hmin : hmax, wmax : wmax) = 255;
        BoxRat = (hmax-hmin)*(wmax-wmin)/(m1*n1);
    end
%     figure; imshow(uint8(Imframe0))
    if NumT<2
        outObjectPath = sprintf('%s/LGFMotionBox%05d.jpg',outPath, j);    
        imwrite(uint8(Imframe0), outObjectPath);
    elseif NumT==2
        outObjectPath = sprintf('%s/FOSMotionBox%05d.jpg',outPath, j);    
        imwrite(uint8(Imframe0), outObjectPath);
    else
        outObjectPath = sprintf('%s/T%dMotionBox%05d.jpg',outPath, NumT, j);    
        imwrite(uint8(Imframe0), outObjectPath);
    end
    
    if Row == 1
        TargetCoordR = [j TargetTempSave BoxRat];
    else   
        TargetNoInside = gzDeleteInsideblock(TargetTempSave);       
        Row = size(TargetNoInside,1);
        for ii = 1:Row
            boxtemp = TargetNoInside(ii,:);
            if ii>1
                boxtempM = repmat(boxtemp(:)', [size(TargetNoOverlap,1) 1]);
                boxtempD = abs(boxtempM - TargetNoOverlap);
                boxtempD = sum(boxtempD,2);
                boxtempDLab = (boxtempD>0);
                if sum(boxtempDLab(:))==(size(TargetNoOverlap,1))
                    TargetNoOverlap0 = [TargetNoOverlap; boxtemp];
                    TargetNoOverlap = gzDeleteInsideblock(TargetNoOverlap0);
                    if size(TargetNoOverlap,1)<size(TargetNoOverlap0,1)
                        continue;
                    end
                end

            else
                TargetNoOverlap = [TargetNoOverlap; boxtemp];
            end

            if ii==1
                for jj = ii+1:Row               
                    boxnext = TargetNoInside(jj,:);
                    box1 = [boxtemp(2),boxtemp(1),(boxtemp(4)-boxtemp(2)),(boxtemp(3)-boxtemp(1))];   %[X,Y,WIDTH,HEIGHT]
                    box2 = [boxnext(2),boxnext(1),(boxnext(4)-boxnext(2)),(boxnext(3)-boxnext(1))];
                    iou = inters_union(box1,box2);
                    if iou > 0
                        minH = min(boxtemp(1),boxnext(1));
                        minW = min(boxtemp(2),boxnext(2));
                        maxH = max(boxtemp(3),boxnext(3));
                        maxW = max(boxtemp(4),boxnext(4));
                        boxcom = [minH minW maxH maxW];
                        TargetNoOverlap(ii,:) = boxcom;
                        boxtemp = boxcom;
                    else
                        TargetNoOverlap = [TargetNoOverlap; boxnext];
                    end               
                end
            else
                Row2 = size(TargetNoOverlap,1);
                for jj = 1:Row2
                    boxnext = TargetNoOverlap(jj,:);
                    box1 = [boxtemp(2),boxtemp(1),(boxtemp(4)-boxtemp(2)),(boxtemp(3)-boxtemp(1))];   %[X,Y,WIDTH,HEIGHT]
                    box2 = [boxnext(2),boxnext(1),(boxnext(4)-boxnext(2)),(boxnext(3)-boxnext(1))];
                    iou = inters_union(box1,box2);
                    if (iou > 0) && (iou < 1)
                        minH = min(boxtemp(1),boxnext(1));
                        minW = min(boxtemp(2),boxnext(2));
                        maxH = max(boxtemp(3),boxnext(3));
                        maxW = max(boxtemp(4),boxnext(4));
                        boxcom = [minH minW maxH maxW];
                        TargetNoOverlap(jj,:) = boxcom;
                    end
                end
                TargetNoOverlap = gzDeleteInsideblock(TargetNoOverlap);
            end     
        end
        TargetNoOverlap = gzDeleteInsideblock(TargetNoOverlap);

        % The second time refinement    
        RowSec = size(TargetNoOverlap,1);
        if RowSec > 2
            TargetNoOverlapSec = [];
            for ii = RowSec:-1:1
                boxtemp = TargetNoOverlap(RowSec,:);
                if ii == 1
                    TargetNoOverlapSec = [TargetNoOverlapSec; TargetNoOverlap(1,:)];
                elseif ii == RowSec
                    TargetNoOverlapSec = [TargetNoOverlapSec; TargetNoOverlap(RowSec,:)];
                else
                    boxnext = TargetNoOverlap(ii,:);
                    box1 = [boxtemp(2),boxtemp(1),(boxtemp(4)-boxtemp(2)),(boxtemp(3)-boxtemp(1))];   %[X,Y,WIDTH,HEIGHT]
                    box2 = [boxnext(2),boxnext(1),(boxnext(4)-boxnext(2)),(boxnext(3)-boxnext(1))];
                    iou = inters_union(box1,box2);
                    if (iou > 0) && (iou < 1)
                        minH = min(boxtemp(1),boxnext(1));
                        minW = min(boxtemp(2),boxnext(2));
                        maxH = max(boxtemp(3),boxnext(3));
                        maxW = max(boxtemp(4),boxnext(4));
                        boxcom = [minH minW maxH maxW];
                        TargetNoOverlapSec = [TargetNoOverlapSec; boxcom];
                    else
                        TargetNoOverlapSec = [TargetNoOverlapSec; boxnext];
                    end
                end
                TargetNoOverlapSec = gzDeleteInsideblock(TargetNoOverlapSec);
            end           
            TargetNoOverlap = TargetNoOverlapSec;
        end     

        % If the detected boxes more than 5, delected the small boxes, only keep 5 boxes
        while size(TargetNoOverlap,1) > 5
            TargetNoOverlapDel = [];
            for iDel = 1:size(TargetNoOverlap,1)
                boxtemp = TargetNoOverlap(iDel,:);
                iArea = (boxtemp(4)-boxtemp(2))*(boxtemp(3)-boxtemp(1));
                TargetNoOverlapDel = [TargetNoOverlapDel;iArea];
            end
            iarMin = find(TargetNoOverlapDel==min(TargetNoOverlapDel));
            TargetNoOverlap(iarMin,:) = [];
        end
        
        if size(TargetNoOverlap,1) <= 3
            % The third time refinement 
            TempImg = zeros(m1, n1);
            for ii = 1:size(TargetNoOverlap,1)
                boxtemp = TargetNoOverlap(ii,:);
                TempImg(boxtemp(1):boxtemp(3),boxtemp(2):boxtemp(4))=1;
            end
            EdgImg = (TempImg>0);
            [L,numtarget]=bwlabel(EdgImg,8);
            TargetNoOverlap = [];
            for nT = 1:numtarget
                [rtemp,ctemp] = find(L==nT);
                minH=min(rtemp);
                minW=min(ctemp);
                maxH=max(rtemp);
                maxW=max(ctemp);

                minH = max(1,minH);
                minW = max(1,minW);
                maxH = min(m1,maxH);
                maxW = min(n1,maxW);
                TargetNoOverlap = [TargetNoOverlap; minH minW maxH maxW];
            end
        end
        
        
        [ms, ns] = size(TargetNoOverlap);
        for iBox = 1:ms,
          blocktemp = TargetNoOverlap(iBox,:);
          TargetCoord = [TargetCoord; j blocktemp];
        end

        Row = size(TargetCoord,1);
        for i = 1:Row
            coortemp = TargetCoord(i,:);
            hmin = coortemp(2);
            wmin = coortemp(3);
            hmax = coortemp(4);
            wmax = coortemp(5);
            BoxRat = (hmax-hmin)*(wmax-wmin)/(m1*n1);
            % Labeling the boxes
            Imframe(hmin : hmin, wmin : wmax) = 255;
            Imframe(hmax : hmax, wmin : wmax) = 255;
            Imframe(hmin : hmax, wmin : wmin) = 255;
            Imframe(hmin : hmax, wmax : wmax) = 255;
            TargetCoordR = [TargetCoordR; coortemp BoxRat];
        end
%         figure; imshow(uint8(Imframe))
        if NumT<2       
            outObjectPath = sprintf('%s/LGFMotionBoxRef%05d.jpg',outPath, j);    
            imwrite(uint8(Imframe), outObjectPath);
        elseif NumT==2
            outObjectPath = sprintf('%s/FOSMotionBoxRef%05d.jpg',outPath, j);    
            imwrite(uint8(Imframe), outObjectPath);
        else
            outObjectPath = sprintf('%s/T%dMotionBoxRef%05d.jpg',outPath, NumT, j);    
            imwrite(uint8(Imframe), outObjectPath);
        end
    end    
    frameObjectBox{j} = TargetCoordR;
    
else
    if j == 1
        TargetCoordR = [j round(0.3*m1) round(0.3*n1) round(0.7*m1) round(0.7*n1) 0];
        frameObjectBox{j} = TargetCoordR;
    else
        frameObjectBox{j} = frameObjectBox{j-1};
    end
    
    coortemp = frameObjectBox{j};
    hmin = coortemp(2);
    wmin = coortemp(3);
    hmax = coortemp(4);
    wmax = coortemp(5);
    % Labeling the boxes
    Imframe(hmin : hmin, wmin : wmax) = 255;
    Imframe(hmax : hmax, wmin : wmax) = 255;
    Imframe(hmin : hmax, wmin : wmin) = 255;
    Imframe(hmin : hmax, wmax : wmax) = 255;
    
%     figure; imshow(uint8(Imframe))
    if NumT<2       
        outObjectPath = sprintf('%s/LGFMotionBoxRef%05d.jpg',outPath, j);    
        imwrite(uint8(Imframe), outObjectPath);
    elseif NumT==2 
        outObjectPath = sprintf('%s/FOSMotionBoxRef%05d.jpg',outPath, j);    
        imwrite(uint8(Imframe), outObjectPath);
    else
        outObjectPath = sprintf('%s/T%dMotionBoxRef%05d.jpg',outPath, NumT, j);    
        imwrite(uint8(Imframe), outObjectPath);
    end
end

if j == nframe-1
  frameObjectBox(j+1) = frameObjectBox(j);
end
ObjectBox = frameObjectBox;


        
%         % The third time refinement 
%         RowThid = size(TargetNoOverlap,1);
%         TargetNoOverlapThid = [];
%         boxtemp = TargetNoOverlap(1,:);
%         for ii = 1:RowThid            
%             if ii == 1
%                 TargetNoOverlapThid = [TargetNoOverlapThid; boxtemp];
%             else
%                 boxnext = TargetNoOverlap(ii,:);
%                 box1 = [boxtemp(2),boxtemp(1),(boxtemp(4)-boxtemp(2)),(boxtemp(3)-boxtemp(1))];   %[X,Y,WIDTH,HEIGHT]
%                 box2 = [boxnext(2),boxnext(1),(boxnext(4)-boxnext(2)),(boxnext(3)-boxnext(1))];
%                 iou = inters_union(box1,box2);
%                 if (iou > 0) && (iou < 1)
%                     minH = min(boxtemp(1),boxnext(1));
%                     minW = min(boxtemp(2),boxnext(2));
%                     maxH = max(boxtemp(3),boxnext(3));
%                     maxW = max(boxtemp(4),boxnext(4));
%                     boxcom = [minH minW maxH maxW];
%                     TargetNoOverlapThid = [TargetNoOverlapThid; boxcom];
%                 else
%                     TargetNoOverlapThid = [TargetNoOverlapThid; boxnext];
%                 end
%             end
%             TargetNoOverlapThid = gzDeleteInsideblock(TargetNoOverlapThid);
%         end
%         TargetNoOverlapThid = gzDeleteInsideblock(TargetNoOverlapThid);
%         TargetNoOverlap = TargetNoOverlapThid; 