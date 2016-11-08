function [ObjectBox, TargetPixelLabels] = TargetBoxLab(framesdata, PixelEasyMask, PixelSeedMask, H, W, outPath)

[m, n] = size(PixelEasyMask);
nH = H;
nW = W;

frameObjectBox = cell(n, 1);
frameTargetPixels = cell(n, 1);

for j = 1:n-1
    TargetCoordR = [];

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
%     index = find(area > 15^2);         % Finding the index of the region that its size larger than a threshold (e.g.10*10)   
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
    index = find(area > 20^2);         % (UCF-->20^2; SegTrack-->5^2)Finding the index of the region that its size larger than a threshold (e.g.10^2)   
    % Obtaining the indexed regions which eliminate small noisy regions 
    TargetBWC = ismember(L,index(:));
    [L,numtarget] = bwlabel(TargetBWC,8);
    TargetPixels = TargetBWC; 
    
    outputori = sprintf('%s/ori%05d.jpg',outPath, j); 
    outputseed  = sprintf('%s/seed%05d.jpg',outPath, j);
    outputTarget  = sprintf('%s/target%05d.jpg',outPath, j);   
    imwrite(orimask, outputori);
    imwrite(seedmask, outputseed);
    imwrite(TargetPixels, outputTarget);    
    
    if numtarget>=1
        frameTargetPixels{j} = TargetBWC;
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
        TargetTempSave = gzDeleteInsideblock(TargetTempSave);
        
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
            TargetCoordR = [TargetCoordR; j coortemp BoxRat];
        end

        outObjectPath = sprintf('%s/OriRPCACleanBoxs%05d.jpg',outPath, j);    
        imwrite(uint8(Imframe0), outObjectPath);

        frameObjectBox{j} = TargetCoordR;

    else
        if j == 1
            TargetCoordR = [j round(0.3*m1) round(0.3*n1) round(0.7*m1) round(0.7*n1) 0];
            frameObjectBox{j} = TargetCoordR;
        else
            frameObjectBox{j}    = frameObjectBox{j-1};
            frameTargetPixels{j} = frameTargetPixels{j-1};
        end
    end
    
end

frameObjectBox(n) = frameObjectBox(n-1);
frameTargetPixels(n) = frameTargetPixels(n-1);

ObjectBox = frameObjectBox;
TargetPixelLabels = frameTargetPixels;

% ObjectsCoord = sprintf('%s/frameTargetBox.mat',outPath);
% save(ObjectsCoord,'ObjectBox');
