function gzFlowNormalizationShow(normflowdata)

[m n] = size(normflowdata);
nW= 4*sqrt(m/12);
nH= 3*sqrt(m/12);

% nSetBlock=nChangeblocksize; 
% nWblock=nW/nSetBlock;
% nHblock=nH/nSetBlock;
for j = 1:n,
    frameTemp   =  normflowdata(:,j);
    maxTemp     =  max(frameTemp);
    frameFlow   = reshape( frameTemp, nH, nW);
    frameNormalize = frameFlow/maxTemp;
    
    outputNormalization  = sprintf('gzNormalizeflow%05d.bmp',j);
    imwrite(frameNormalize, outputNormalization);
end 

end
