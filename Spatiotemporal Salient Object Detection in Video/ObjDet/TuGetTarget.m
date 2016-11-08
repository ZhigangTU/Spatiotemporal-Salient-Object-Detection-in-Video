function TargetCoord = TuGetTarget(PixelEasyMask, PixelSeedMask)

[m, n] = size(PixelEasyMask);
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);

% TargetCoord = zeros( n , 4 );
TargetCoord = [];
TargetCoordPredict = [];

for j = 1:n,
  orimask  = reshape(PixelEasyMask(:,j), nH, nW);
  seedmask = reshape(PixelSeedMask(:,j), nH, nW);
  
  outputori = sprintf('zori%05d.bmp',j);
  imwrite(orimask, outputori);
  outputseed  = sprintf('zseed%05d.bmp',j);
  imwrite(seedmask, outputseed);
  
  %to grow by using the pixel seed
  TargetPixels = gzSeedgrow(seedmask,orimask);
  outputTarget = sprintf('ztarget%05d.bmp',j);
  imwrite(TargetPixels, outputTarget);
  
  % try to use some matlab functions
  TargetBW = TargetPixels>0.6;
  [L,numtarget]=bwlabel(TargetBW,8);  % original(8); 4
  [m1, n1] = size(TargetPixels);
  
  if numtarget>=1
      TargetCoordPredict = [];
      TargetTempSave = [];
      
      if numtarget == 1
          [rtemp,ctemp] = find(L==1);
          minH=min(rtemp);
          minW=min(ctemp);
          maxH=max(rtemp);
          maxW=max(ctemp);
          
          nExtend=1; % nExtend=2;
          minH = max(1,minH - nExtend);
          minW = max(1,minW - nExtend);
          maxH = min(m1,maxH + nExtend);
          maxW = min(n1,maxW + nExtend);
          TargetTempSave = [TargetTempSave; minH minW maxH maxW];
          TargetCoord = [TargetCoord; j minH minW maxH maxW];
          TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW];
          
      elseif numtarget == 2
          [rtemp1,ctemp1] = find(L==1);
          [rtemp2,ctemp2] = find(L==2);
          B1minH=min(rtemp1);
          B1minW=min(ctemp1);
          B1maxH=max(rtemp1);
          B1maxW=max(ctemp1);
          B2minH=min(rtemp2);
          B2minW=min(ctemp2);
          B2maxH=max(rtemp2);
          B2maxW=max(ctemp2);
          
          I1 = zeros(nH, nW);
          I2 = zeros(nH, nW);
          I1(B1minH:B1maxH, B1minW:B1maxW) = 1;
          I2(B2minH:B2maxH, B2minW:B2maxW) = 1;
          II = I1.*I2;
          if sum(II(:))>=1
              minH = min(B1minH,B2minH);
              minW = min(B1minW,B2minW);
              maxH = max(B1maxH,B2maxH);
              maxW = max(B1maxW,B2maxW);
              TargetTempSave = [TargetTempSave; minH minW maxH maxW];
              TargetCoord = [TargetCoord; j minH minW maxH maxW];
              TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW];
          else
              TargetCoord = [TargetCoord; j B1minH B1minW B1maxH B1maxW];
              TargetCoord = [TargetCoord; j B2minH B2minW B2maxH B2maxW];
              TargetCoordPredict = [TargetCoordPredict; j+1 B2minH B2minW B2maxH B2maxW];
          end
          
      elseif numtarget == 3
          [rtemp1,ctemp1] = find(L==1);
          [rtemp2,ctemp2] = find(L==2);
          [rtemp3,ctemp3] = find(L==3);
          B1minH=min(rtemp1);
          B1minW=min(ctemp1);
          B1maxH=max(rtemp1);
          B1maxW=max(ctemp1);
          B2minH=min(rtemp2);
          B2minW=min(ctemp2);
          B2maxH=max(rtemp2);
          B2maxW=max(ctemp2);
          B3minH=min(rtemp3);
          B3minW=min(ctemp3);
          B3maxH=max(rtemp3);
          B3maxW=max(ctemp3);
          I1 = zeros(nH, nW);
          I2 = zeros(nH, nW);
          I3 = zeros(nH, nW);
          I1(B1minH:B1maxH, B1minW:B1maxW) = 1;
          I2(B2minH:B2maxH, B2minW:B2maxW) = 1;
          I3(B3minH:B3maxH, B3minW:B3maxW) = 1;
          II12 = I1.*I2;
          II13 = I1.*I3;
          II23 = I2.*I3;
          
          if sum(II12(:))>=1
              minH = min(B1minH,B2minH);
              minW = min(B1minW,B2minW);
              maxH = max(B1maxH,B2maxH);
              maxW = max(B1maxW,B2maxW);
              IIC12 = zeros(nH, nW);
              IIC12(minH:maxH, minW:maxW) = 1;
              
              IIC123 = IIC12.*I3;
              if sum(IIC123(:))>=1
                  minH = min(minH,B3minH);
                  minW = min(minW,B3minW);
                  maxH = max(maxH,B3maxH);
                  maxW = max(maxW,B3maxW);
                  TargetTempSave = [TargetTempSave; minH minW maxH maxW];
                  TargetCoord = [TargetCoord; j minH minW maxH maxW];
                  TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW];
              else
                  TargetCoord = [TargetCoord; j minH minW maxH maxW];
                  TargetCoord = [TargetCoord; j B3minH B3minW B3maxH B3maxW];
                  TargetCoordPredict = [TargetCoordPredict; j+1 B3minH B3minW B3maxH B3maxW];
              end
          elseif sum(II13(:))>=1
              minH = min(B1minH,B3minH);
              minW = min(B1minW,B3minW);
              maxH = max(B1maxH,B3maxH);
              maxW = max(B1maxW,B3maxW);
              IIC13 = zeros(nH, nW);
              IIC13(minH:maxH, minW:maxW) = 1;
                  
              IIC123 = IIC13.*I2;
              if sum(IIC123(:))>=1
                  minH = min(minH,B2minH);
                  minW = min(minW,B2minW);
                  maxH = max(maxH,B2maxH);
                  maxW = max(maxW,B2maxW);
                  TargetTempSave = [TargetTempSave; minH minW maxH maxW];
                  TargetCoord = [TargetCoord; j minH minW maxH maxW];
                  TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW];
              else
                  TargetCoord = [TargetCoord; j minH minW maxH maxW];
                  TargetCoord = [TargetCoord; j B2minH B2minW B2maxH B2maxW];
                  TargetCoordPredict = [TargetCoordPredict; j+1 B2minH B2minW B2maxH B2maxW];
              end
          elseif sum(II23(:))>=1
              minH = min(B2minH,B3minH);
              minW = min(B2minW,B3minW);
              maxH = max(B2maxH,B3maxH);
              maxW = max(B2maxW,B3maxW);
              IIC23 = zeros(nH, nW);
              IIC23(minH:maxH, minW:maxW) = 1;
              
              IIC123 = IIC23.*I1;
              if sum(IIC123(:))>=1
                  minH = min(minH,B1minH);
                  minW = min(minW,B1minW);
                  maxH = max(maxH,B1maxH);
                  maxW = max(maxW,B1maxW);
                  TargetTempSave = [TargetTempSave; minH minW maxH maxW];
                  TargetCoord = [TargetCoord; j minH minW maxH maxW];
                  TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW];
              else
                  TargetCoord = [TargetCoord; j minH minW maxH maxW];
                  TargetCoord = [TargetCoord; j B1minH B1minW B1maxH B1maxW];
                  TargetCoordPredict = [TargetCoordPredict; j+1 B1minH B1minW B1maxH B1maxW];
              end
              
          else
              TargetCoord = [TargetCoord; j B1minH B1minW B1maxH B1maxW];
              TargetCoord = [TargetCoord; j B2minH B2minW B2maxH B2maxW];
              TargetCoord = [TargetCoord; j B3minH B3minW B3maxH B3maxW];
              TargetCoordPredict = [TargetCoordPredict; j+1 B3minH B3minW B3maxH B3maxW];
          end
      end
  else
      TargetCoordPredict(1)=j;
      TargetCoord = [TargetCoord;TargetCoordPredict];
  end
                                    
%       if numtarget>=1
%           TargetCoordPredict = [];
%           TargetTempSave = [];
%           for nT = 1:numtarget,
%               [rtemp,ctemp] = find(L==nT);
%               minH=min(rtemp);
%               minW=min(ctemp);
%               maxH=max(rtemp);
%               maxW=max(ctemp);
% 
%               %nExtend=2;
%               nExtend=1;
%               minH = max(1,minH - nExtend);
%               minW = max(1,minW - nExtend);
%               maxH = min(m1,maxH + nExtend);
%               maxW = min(n1,maxW + nExtend);
%               TargetTempSave = [TargetTempSave; minH minW maxH maxW];
%               TargetCoord = [TargetCoord; j minH minW maxH maxW];
%               TargetCoordPredict = [TargetCoordPredict; j+1 minH minW maxH maxW ];
%           end
% %           TargetNoInside = gzDeleteInsideblock(TargetTempSave);
% %           [ms ns] = size(TargetNoInside);
% %           for ngood = 1:ms,
% %               blocktemp = TargetNoInside(ngood,:);
% %               TargetCoord = [TargetCoord; j blocktemp];
% %               TargetCoordPredict = [TargetCoordPredict; j+1 blocktemp];
% %           end
%       else
%           TargetCoord = [TargetCoord;TargetCoordPredict];
%       end
  
  
%   %treeman case, i write it by myself function
%   [minH minW maxH maxW]  = gzGetCoordofTarget(TargetPixels);
%   TargetCoord(j,:) = [minH minW maxH maxW];
%   if (minH>maxH|minW>maxW)&(j>1)
%       TargetCoord(j,:) = TargetCoord(j-1,:)
%   end
 end

end
