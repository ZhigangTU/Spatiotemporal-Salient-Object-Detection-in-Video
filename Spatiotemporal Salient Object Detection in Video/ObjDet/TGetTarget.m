function [TargetCoord, TargetMask]= TGetTarget(PixelEasyMask, PixelSeedMask, H, W, outPath)

[m, n] = size(PixelEasyMask);
TargetMask = zeros(m, n);
% nW= 4*sqrt(m/12);
% nH= 3*sqrt(m/12);
nH = H;
nW = W;

% TargetCoord = zeros(n , 4);
TargetCoord = [];
TargetCoordPredict = [];
for j = 1:n    
    orimask  = reshape(PixelEasyMask(:,j), nH, nW);
    seedmask = reshape(PixelSeedMask(:,j), nH, nW);
  
    % To grow by using the pixel seed
    if sum(seedmask(:))>1
        TargetPixels = gzSeedgrow(seedmask,orimask);  
    else
        TargetPixels = orimask;
    end
  
    % Try to use some matlab functions
    TargetBW = TargetPixels>0.6;
    TargetBW = EdgeClean(TargetBW, nH, nW, 1);
    [L,numtarget]=bwlabel(TargetBW,8);
    [m1, n1] = size(TargetPixels);

    % Delect and clean noisy pixels
    stats = regionprops(L,'Area');  % Compute the size of each connected region
    area = cat(1,stats.Area);     
    index = find(area > 100);       % Finding the index of the region that its size larger than a threshold (e.g.10*10)   
    % Obtaining the indexed regions which eliminate small noisy regions 
    TargetBWC = ismember(L,index(:));
    [L,numtarget]=bwlabel(TargetBWC,8);
    TargetPixels = TargetBWC; 
  
   % Record the objects Labels
    TargetMask(:,j) = reshape(TargetPixels, m, 1);
  
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
          
            % nExtend=2;
            nExtend=1;
            minH = max(1,minH - nExtend);
            minW = max(1,minW - nExtend);
            maxH = min(m1,maxH + nExtend);
            maxW = min(n1,maxW + nExtend);
            TargetTempSave = [TargetTempSave; minH minW maxH maxW];
            % TargetCoord = [TargetCoord; j minH minW maxH maxW];
            % TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW ];
        end
        TargetNoInside = gzDeleteInsideblock(TargetTempSave);
        [ms, ns] = size(TargetNoInside);
        for ngood = 1:ms,
            blocktemp = TargetNoInside(ngood,:);
            TargetCoord = [TargetCoord; j blocktemp];
            TargetCoordPredict = [TargetCoordPredict; j+1 blocktemp];
        end
    else
        TargetCoord = [TargetCoord;TargetCoordPredict];
    end
  
    if j == n
        TargetCoord(j) = TargetCoord(j-1);
    end
   
end

%   %treeman case, i write it by myself function
%   [minH minW maxH maxW]  = gzGetCoordofTarget(TargetPixels);
%   TargetCoord(j,:) = [minH minW maxH maxW];
%   if (minH>maxH|minW>maxW)&(j>1)
%       TargetCoord(j,:) = TargetCoord(j-1,:)
%   end
