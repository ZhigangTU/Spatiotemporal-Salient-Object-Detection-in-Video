function SelectBlocks(data, TargetCoord)

[m,n] = size(data);
[m1,n1] = size(TargetCoord);
nW = 4*sqrt(m/12);
nH = 3*sqrt(m/12);
for j = 1:n
  imgframe = reshape(data(:,j), nH, nW);
  
  for k = 1:m1
      coortemp = TargetCoord(k,:);
      if coortemp(1) == j
           hmin = coortemp(2);
           wmin = coortemp(3);
           hmax = coortemp(4);
           wmax = coortemp(5);
           
           
           
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
  
  outputNum  = sprintf('ztagetblcok%05d.bmp',j);
  outputimg  = fullfile (destDir, outputNum);
  imwrite(imgframe, outputimg);
end



%to grow by using the pixel seed
           imgframe(hmin : hmin, wmin : wmax) = 1;
           imgframe(hmax : hmax, wmin : wmax) = 1;
           imgframe(hmin : hmax, wmin : wmin) = 1;
           imgframe(hmin : hmax, wmax : wmax) = 1;









