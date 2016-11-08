function ObjectBox = TargetBoxLab(framesdata, PixelEasyMask, PixelSeedMask, H, W, outPath, NumT)

if nargin < 6
    NumT = 1;
end

[m, n] = size(PixelEasyMask);
nH = H;
nW = W;

frameObjectBox = cell(n, 1);

for j = 1:n-1
    TargetCoord = [];
    TargetCoordR = [];
    TargetNoOverlap = [];

    Imframe = framesdata{j};
    
    orimask  = reshape(PixelEasyMask(:,j), nH, nW);
    seedmask = reshape(PixelSeedMask(:,j), nH, nW);
  
    % To grow by using the pixel seed
    if sum(seedmask(:))>1
        TargetPixels = gzSeedgrow(seedmask,orimask);  
    else
        TargetPixels = orimask;
    end
  
%     % % 1) JHMDB
%     nExtend = 2;     % JHMDB-->2; Sports-->3
%     [m1, n1] = size(TargetPixels);
%     TargetBW = TargetPixels>0.6;
%     % Clean noisy pixels
%     TargetBW = EdgeClean(TargetBW, nH, nW, 1); 
%     [L,numtarget]=bwlabel(TargetBW,8);
%     
%     % Delect and clean noisy pixels
%     stats = regionprops(L,'Area');    % Compute the size of each connected region
%     area = cat(1,stats.Area);     
%     index = find(area > 100);         % Finding the index of the region that its size larger than a threshold (e.g.10*10)   
%     % Obtaining the indexed regions which eliminate small noisy regions 
%     TargetBWC = ismember(L,index(:));
%     [L,numtarget]=bwlabel(TargetBWC,8);
%     TargetPixels = TargetBWC; 
    
    % % 2) UCF-Sports
    nExtend = 3;     % JHMDB-->2; Sports-->3
    [m1, n1] = size(TargetPixels);
    TargetBW = TargetPixels>0.6;
    % Clean noisy pixels
    TargetBW = EdgeClean(TargetBW, nH, nW, 1);
    % Dilate the Labeled objects
    TargetBWD = imdilate(TargetBW,strel('diamond',1));    
    [L,numtarget]=bwlabel(TargetBWD,8);
  
    % Delect and clean noisy pixels
    stats = regionprops(L,'Area');     % Compute the size of each connected region
    area = cat(1,stats.Area);     
    index = find(area > 20^2);         % Finding the index of the region that its size larger than a threshold (e.g.10^2)   
    % Obtaining the indexed regions which eliminate small noisy regions 
    TargetBWC = ismember(L,index(:));
    [L,numtarget] = bwlabel(TargetBWC,8);
    TargetPixels = TargetBWC; 
    
    outputori = sprintf('%s/ori%05d.bmp',outPath, j); 
    outputseed  = sprintf('%s/seed%05d.bmp',outPath, j);
    outputTarget  = sprintf('%s/target%05d.bmp',outPath, j);   
    imwrite(orimask, outputori);
    imwrite(seedmask, outputseed);
    imwrite(TargetPixels, outputTarget); 
    
    if numtarget>=1
        TargetTempSave = [];
        for nT = 1:numtarget,
            [rtemp,ctemp] = find(L==nT);
            minH=min(rtemp);
            minW=min(ctemp);
            maxH=max(rtemp);
            maxW=max(ctemp);

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
%         figure; imshow(uint8(Imframe0))
        if NumT<2
            outObjectPath = sprintf('%s/ObjectBox%05d.bmp',outPath, j);    
            imwrite(uint8(Imframe0), outObjectPath);
        else
            outObjectPath = sprintf('%s/T%dObjectBox%05d.bmp',outPath, NumT, j);    
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
            
%             Row = size(TargetNoOverlap,1);
%             for i = 1:Row
%                 coortemp = TargetNoOverlap(i,:);
%                 hmin = coortemp(1);
%                 wmin = coortemp(2);
%                 hmax = coortemp(3);
%                 wmax = coortemp(4);
%                 % Labeling the boxes
%                 Imframe(hmin : hmin, wmin : wmax) = 255;
%                 Imframe(hmax : hmax, wmin : wmax) = 255;
%                 Imframe(hmin : hmax, wmin : wmin) = 255;
%                 Imframe(hmin : hmax, wmax : wmax) = 255;
%             end
%             figure; imshow(uint8(Imframe))           

            % The second time refinement    
            RowSec = size(TargetNoOverlap,1);
            if RowSec > 2
                TargetNoOverlapSec = [];
                boxtemp = TargetNoOverlap(RowSec,:);
                for ii = RowSec:-1:1                    
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
            
            % The third time refinement    
            RowThid = size(TargetNoOverlap,1);
            if RowThid > 3
                TargetNoOverlapThid = [];
                iRow = 2;
                boxtemp = TargetNoOverlap(iRow,:);
                TargetNoOverlapThid = [TargetNoOverlapThid; TargetNoOverlap(iRow,:)];
                for ii = 1:RowThid                 
                    if ii == iRow   
                        continue;
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
                            TargetNoOverlapThid = [TargetNoOverlapThid; boxcom];
                        else
                            TargetNoOverlapThid = [TargetNoOverlapThid; boxnext];
                        end
                    end
                    TargetNoOverlapThid = gzDeleteInsideblock(TargetNoOverlapThid);
                end
                TargetNoOverlapThid = gzDeleteInsideblock(TargetNoOverlapThid);
                TargetNoOverlap = TargetNoOverlapThid;
            end

            [ms, ns] = size(TargetNoOverlap);
            for iBox = 1:ms
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
%             figure; imshow(uint8(Imframe))
            if NumT<2       
                outObjectPath = sprintf('%s/ObjectBoxRef%05d.bmp',outPath, j);    
                imwrite(uint8(Imframe), outObjectPath);
            else
                outObjectPath = sprintf('%s/T%dObjectBoxRef%05d.bmp',outPath, NumT, j);    
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
    end    
end

frameObjectBox(n) = frameObjectBox(n-1);

ObjectBox = frameObjectBox;
    