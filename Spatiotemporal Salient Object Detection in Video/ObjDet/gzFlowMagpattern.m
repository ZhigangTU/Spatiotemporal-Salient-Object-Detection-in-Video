function FlowPattern = gzFlowMagpattern(normflowdata, nChangeblocksize)

[m n] = size(normflowdata);
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);

nSetBlock=nChangeblocksize; 
%nSetBlock=4; %the block size will be nSetBlock*nSetBlock
nWblock=nW/nSetBlock;
nHblock=nH/nSetBlock;

FlowPattern = zeros( nWblock*nHblock, n);

for j = 1:n,
  PatternTemp = zeros( nHblock, nWblock);
  frameTemp   = reshape( normflowdata(:,j),nH,nW);
  
  for jh = 2:nHblock-1,
  for jw = 2:nWblock-1,
      %flow norm of central block
      nhStart=(jh-1)*nSetBlock+1;
      nwStart=(jw-1)*nSetBlock+1;
      BlockTemp = frameTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1);
      Blocknorm = norm(BlockTemp,'fro');
      
      xNeibor=[1 1 0 -1 -1 -1 0 1];
      yNeibor=[0 1 1 1 0 -1 -1 -1];
      NeiborSum = 0;
      for nAdd = 1:8,
          neiH = jh + yNeibor(nAdd);
          neiW = jw + xNeibor(nAdd);
          
          nhNStart=(neiH-1)*nSetBlock+1;
          nwNStart=(neiW-1)*nSetBlock+1;
          
          NeiborTemp = frameTemp(nhNStart : nhNStart + nSetBlock -1, nwNStart : nwNStart + nSetBlock -1);
          Neibornorm = norm(NeiborTemp,'fro');
          NeiborSum = NeiborSum + Neibornorm;
      end
      NeiborAve = NeiborSum/8;
      
      ratio = Blocknorm/NeiborAve;
      Pvalue =0;
%       if ratio<=0.7
%           Pvalue = 1;
%       elseif (ratio>0.7) & (ratio<1.3)
%           Pvalue = 0;%0.5;
%       elseif ratio>=1.3
%           Pvalue = 1;
%       end 
      
      if ratio>=0.85
          Pvalue = 1;
      end 
      PatternTemp(jh,jw) = Pvalue; 
  end
  end
  
  Patternvector = reshape(PatternTemp,nHblock*nWblock,1);
  FlowPattern(:,j)=Patternvector;
  
  ShowTemp  = zeros(nH,nW);
  for jh = 1:nHblock,
  for jw = 1:nWblock,
      nhStart=(jh-1)*nSetBlock+1;
      nwStart=(jw-1)*nSetBlock+1;
      ShowTemp(nhStart : nhStart + nSetBlock -1, nwStart : nwStart + nSetBlock -1)=PatternTemp(jh,jw)*ones(nSetBlock,nSetBlock);
  end
  end
  
  outputPattern  = sprintf('gzflowmag%05d.bmp',j);
  imwrite(ShowTemp, outputPattern);
end 


end
