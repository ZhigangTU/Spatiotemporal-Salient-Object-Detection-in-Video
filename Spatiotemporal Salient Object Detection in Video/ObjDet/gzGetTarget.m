function [TargetCoord, TargetMask]= gzGetTarget(PixelEasyMask, PixelSeedMask, H, W, outPath, )

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
  
  outputori = sprintf('%s/ori%05d.bmp',outPath, j); 
  imwrite(orimask, outputori);
  outputseed  = sprintf('%s/seed%05d.bmp',outPath, j);
  imwrite(seedmask, outputseed);
  
  %to grow by using the pixel seed
  TargetPixels = gzSeedgrow(seedmask,orimask);
  
  TargetMask(:,j) = reshape(TargetPixels, m, 1);
  
  outputTarget  = sprintf('%s/target%05d.bmp',outPath, j);
  imwrite(TargetPixels, outputTarget);
  
  
  
  % try to use some matlab functions
  TargetBW = TargetPixels>0.6;
  [L,numtarget]=bwlabel(TargetBW,8);
  [m1 n1] = size(TargetPixels);
  
  
  
  if numtarget>=1
      TargetCoordPredict = [];
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
      [ms ns] = size(TargetNoInside);
      for ngood = 1:ms,
          blocktemp = TargetNoInside(ngood,:);
          TargetCoord = [TargetCoord; j blocktemp];
          TargetCoordPredict = [TargetCoordPredict; j+1 blocktemp];
      end
  else
      TargetCoord = [TargetCoord;TargetCoordPredict];
  end
  
  
%   %treeman case, i write it by myself function
%   [minH minW maxH maxW]  = gzGetCoordofTarget(TargetPixels);
%   TargetCoord(j,:) = [minH minW maxH maxW];
%   if (minH>maxH|minW>maxW)&(j>1)
%       TargetCoord(j,:) = TargetCoord(j-1,:)
%   end
end

end