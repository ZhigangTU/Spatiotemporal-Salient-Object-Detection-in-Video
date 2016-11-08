function PixelMaskRelativeChange = gzComputeRelativeChangePixelmask(data, A_hat, dRatio)
%[A_hat E_hat iter] = inexact_alm_rpca(D, lambda, tol, maxIter)

[m n] = size(data);
RelativeChange = zeros( m, n);

Change         =  abs(data - A_hat);
RelativeChange =  Change./data;

%%
D     = RelativeChange;
%position of each pixel is set with a mask value 1 or 0
PixelMaskRelativeChange   = zeros(m, n);

%gaozhi add
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);


for j = 1:n,
  frameTemp1= reshape( D(:,j),nH,nW);
   
  E_vector   = reshape(frameTemp1,nH*nW,1);
  Esortpixel = sort(E_vector,'descend');
  dThreshpixel  =  Esortpixel(nH*nW*dRatio);%-0.001; %dThresh  = Esort(nHblock*nWblock*0.05);
  if dThreshpixel<0.00001
      dThreshpixel=0.00001;
  end 
  
  frameTemp1    =  frameTemp1>dThreshpixel;
  
  outputMask  = sprintf('zRCpixelmask%05d.bmp',j);
  imwrite(frameTemp1, outputMask);
  PixelMaskRelativeChange(:,j)= reshape(frameTemp1,nH*nW,1);
  
end
%%

end