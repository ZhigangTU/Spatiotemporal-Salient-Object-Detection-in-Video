function [MaskBlock1Pixle2,E_blockmask,PixelMask]=gzEmatrixFindMasksChangingblocksize(NoiseMatrix, dRatio,nChangeblocksize)
%[MaskBlock1Pixle2 E_blockmask]  = gzAnalyzeEmatrixFindMasksChangingblocksize(E_hat,0.3);

currentPath = cd;
tempName = 'mask all kinds of masks';
destDir = fullfile (currentPath,tempName) ;
if ~exist(destDir,'dir')
    mkdir (currentPath,tempName) ;
end

D     = NoiseMatrix;
[m,n] = size(D);

% position of each pixel is set with a mask value 1 or 0
Mask_matrix = zeros(m, n); 
PixelMask   = zeros(m, n);

nW = 4*sqrt(m/12);
nH = 3*sqrt(m/12);
nSetBlock = nChangeblocksize; 
% nSetBlock=4; %the block size will be nSetBlock*nSetBlock
nWblock=nW/nSetBlock;
nHblock=nH/nSetBlock;

% position of each block is set with a mask value 1 or 0
E_blockmask=zeros( nWblock*nHblock, n);
MaskBlock1Pixle2=zeros(m, n);

for j = 1:n,
  E_energy = zeros( nHblock, nWblock);
  
  maskTemp  = reshape( Mask_matrix(:,j),nH,nW);
  frameTemp = reshape( D(:,j),nH,nW);
  frameTemp1= reshape( D(:,j),nH,nW);
  for jh = 1:nHblock,
  for jw = 1:nWblock,
      nhStart=(jh-1)*nSetBlock+1;
      nwStart=(jw-1)*nSetBlock+1;
      BlockTemp = frameTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1);
      Blocknorm = norm(BlockTemp,'fro');
      E_energy(jh,jw) = Blocknorm;
  end
  end
  E_energyvector = reshape(E_energy,nHblock*nWblock,1);
  Esort    = sort(E_energyvector,'descend');
  dThresh  = Esort(nHblock*nWblock*dRatio);%-0.001; %dThresh  = Esort(nHblock*nWblock*0.05);
  if dThresh<0.00001
      dThresh=0.00001;
  end 
  
  E_energy    = E_energy>dThresh;
  Eblock      = reshape(E_energy,nHblock*nWblock,1);
  E_blockmask(:,j)=Eblock;
  
  for jh = 1:nHblock
  for jw = 1:nWblock
      nhStart=(jh-1)*nSetBlock+1;
      nwStart=(jw-1)*nSetBlock+1;
      maskTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1)=E_energy(jh,jw)*ones(nSetBlock,nSetBlock);
  end
  end
  
%   outputMasktemp = sprintf('gzblockmask%05d.bmp',j);
%   outputMask  = fullfile (destDir, outputMasktemp);
%   imwrite(maskTemp, outputMask);
  
  % pixel mask
  frameTemp1 = abs(frameTemp1);  
  E_vector   = reshape(frameTemp1,nH*nW,1);
  Esortpixel = sort(E_vector,'descend');
  dThreshpixel  =  Esortpixel(nH*nW*dRatio);%-0.001; %dThresh  = Esort(nHblock*nWblock*0.05);
  if dThreshpixel<0.00001
      dThreshpixel=0.00001;
  end 
  
  frameTemp1 = frameTemp1>dThreshpixel;
%   outputMasktemp= sprintf('gzpixelmask%05d.bmp',j);
%   outputMask    = fullfile (destDir,outputMasktemp);
%   imwrite(frameTemp1, outputMask);
  PixelMask(:,j) = reshape(frameTemp1,nH*nW,1);
  
  MaskAnd       =maskTemp&frameTemp1;
  MaskBlockPixel=(MaskAnd+maskTemp)/2;
%   outputMasktemp= sprintf('BlockPixelmask%05d.bmp',j);
%   outputMask    = fullfile (destDir,outputMasktemp);
%   imwrite(MaskBlockPixel, outputMask);
  
  Mask_temp            = reshape(MaskBlockPixel,nH*nW,1);
  MaskBlock1Pixle2(:,j)= Mask_temp;
end
