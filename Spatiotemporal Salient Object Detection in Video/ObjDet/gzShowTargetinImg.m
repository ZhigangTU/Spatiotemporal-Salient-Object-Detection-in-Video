function gzShowTargetinImg(data, TargetCoord, userName, H, W, Tim)

if nargin < 6
    Tim = 1;
end
        
currentPath = cd;
tempName = 'Target with boxes';
destDir = fullfile (currentPath,tempName, userName) ;
if ~exist(destDir,'dir')
    mkdir (destDir);
end

[m, n] = size(data);
[m1,n1] = size(TargetCoord);
nH = H;
nW = W;

for j = 1:n,
  imgframe  = reshape( data(:,j), nH, nW );
  for k = 1:m1,
      coortemp = TargetCoord(k,:);
      if coortemp(1)==j
           hmin = coortemp(2);
           wmin = coortemp(3);
           hmax = coortemp(4);
           wmax = coortemp(5);
           % to grow by using the pixel seed
           imgframe(hmin : hmin, wmin : wmax) = 1;
           imgframe(hmax : hmax, wmin : wmax) = 1;
           imgframe(hmin : hmax, wmin : wmin) = 1;
           imgframe(hmin : hmax, wmax : wmax) = 1;
      end 
  end
  
%   coortemp  = TargetCoord(j,:);
%   hmin = coortemp(1);
%   wmin = coortemp(2);
%   hmax = coortemp(3);
%   wmax = coortemp(4);
%   %to grow by using the pixel seed
%   imgframe(hmin : hmin, wmin : wmax) = 1;
%   imgframe(hmax : hmax, wmin : wmax) = 1;
%   imgframe(hmin : hmax, wmin : wmin) = 1;
%   imgframe(hmin : hmax, wmax : wmax) = 1;
  if nargin < 6
    outputNum = sprintf('Tagetboxes%05d.bmp', j);
  else
    outputNum = sprintf('%dTagetboxes%05d.bmp', Tim, j);
  end
  outputimg = fullfile (destDir, outputNum);
  imwrite(imgframe, outputimg);
end
