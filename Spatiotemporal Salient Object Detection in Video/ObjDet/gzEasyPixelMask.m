function PixelEasyMask = gzEasyPixelMask(E_hat, dRatio)

currentPath = cd;
tempName = 'mask easy pixel masks';
destDir = fullfile (currentPath,tempName) ;
if ~exist(destDir,'dir')
    mkdir (currentPath,tempName) ;
end

%%
[m, n] = size(E_hat);
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);
%%
RelativeChange =  abs(E_hat);
D     = RelativeChange;
%position of each pixel is set with a mask value 1 or 0
PixelEasyMask = zeros( m, n);

for j = 1:n,
  frameTemp1= reshape( D(:,j),nH,nW);
  E_vector   = reshape(frameTemp1,nH*nW,1);
  Esortpixel = sort(E_vector,'descend');
  dThreshpixel = Esortpixel(nH*nW*dRatio);%-0.001; %dThresh  = Esort(nHblock*nWblock*0.05);
  if dThreshpixel < 0.00001
      dThreshpixel = 0.00001;
  end 
  
  frameTemp1     = frameTemp1 > dThreshpixel;
  outputMasktemp = sprintf('zEasypixelmask%05d.bmp',j);
  outputMask     = fullfile (destDir, outputMasktemp);
  imwrite(frameTemp1, outputMask);
  PixelEasyMask(:,j)= reshape(frameTemp1,nH*nW,1);
end
%%

end