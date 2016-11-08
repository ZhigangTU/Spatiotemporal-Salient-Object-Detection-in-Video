function [Mask_matrix E_blockmask] = gzAnalyzeEmatrix(NoiseMatrix, H, W)
% D - m x n matrix of observations/data (required input)
%
addpath PROPACK;

D = NoiseMatrix;
[m n] = size(D);

Mask_matrix = zeros(m, n);

% %gaozhi add
% nW= 4*sqrt(m/12);
% nH= 3*sqrt(m/12);

nSetBlock=8;       %the block size will be nSetBlock*nSetBlock
nWblock=nW/nSetBlock;
nHblock=nH/nSetBlock;

E_blockmask=zeros( nWblock*nHblock, n);

for j = 1:n,
  E_energy = zeros( nHblock, nWblock);
  
  maskTemp  = reshape( Mask_matrix(:,j),nH,nW);
  frameTemp = reshape( D(:,j),nH,nW);
  for jh = 1:nHblock,
  for jw = 1:nWblock,
      nhStart=(jh-1)*nSetBlock+1;
      nwStart=(jw-1)*nSetBlock+1;
      BlockTemp = frameTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1);
      Blocknorm = norm(BlockTemp,'fro');%/nSetBlock;
      E_energy(jh,jw) = Blocknorm;
  end
  end
  E_energyvector = reshape(E_energy,nHblock*nWblock,1);
  Esort    = sort(E_energyvector,'descend');
  dThresh  = Esort(nHblock*nWblock*1);
  %dThresh  = Esort(nHblock*nWblock*0.05);
  E_energy    = E_energy>=dThresh;
  Eblock      = reshape(E_energy,nHblock*nWblock,1);
  E_blockmask(:,j)=Eblock;
  
  % further process the 8*8 block into 8*8 sub-blocks
  E_sub8by8energy=zeros( nHblock*nSetBlock, nWblock*nSetBlock);
  for jh = 1:nHblock,
  for jw = 1:nWblock,
      if E_energy(jh,jw) == 1
          nhStart=(jh-1)*nSetBlock+1;
          nwStart=(jw-1)*nSetBlock+1;
          OriNoiseTemp = frameTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1);
          OriNoiseTemp = abs(OriNoiseTemp);
          E_sub8by8energy(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1)=OriNoiseTemp;
      end
  end
  end
  E_subvector  = reshape(E_sub8by8energy,nHblock*nWblock*nSetBlock*nSetBlock,1);
  Esubsort     = sort(E_subvector,'descend');
  dSubThresh   = Esubsort(nHblock*nWblock*nSetBlock*nSetBlock*1*0.1);
  Esubmask1    = E_sub8by8energy>0;
  Esubmask2    = E_sub8by8energy>dSubThresh;
  Esubmask     = Esubmask1+Esubmask2;
  
  Esubimage    =  Esubmask * 0.5;
  %imwrite(Esubimage, 'gzsubmask.jpg');
  %imwrite(Esubimage, 'gzsubmask.bmp');
  
  outputMask  = sprintf('gzsubmask%05d.bmp',j);
  imwrite(Esubimage, outputMask);
  %outputMasks = fullfile(destDir, outputMask); 
  
  for jh = 1:nHblock,
  for jw = 1:nWblock,
      nhStart=(jh-1)*nSetBlock+1;
      nwStart=(jw-1)*nSetBlock+1;
      maskTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1)=E_energy(jh,jw)*ones(nSetBlock,nSetBlock);
  end
  end
  maskTemp = reshape(maskTemp,nW*nH,1);
  Mask_matrix(:,j) = maskTemp;
end

%if mod( total_svd, 10) == 0
%        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
%            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
%            ' stopCriterion ' num2str(stopCriterion)]);
%end   
