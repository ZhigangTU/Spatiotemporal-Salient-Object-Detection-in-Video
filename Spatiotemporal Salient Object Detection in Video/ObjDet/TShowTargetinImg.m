function TShowTargetinImg(frameObject, outPath, H, W)

outPath = fullfile(outPath,'RPCA');
if ~exist(outPath,'dir')
    mkdir (outPath);
end

[m, n] = size(frameObject);
nH = H;
nW = W;

for j = 1:m
    ObjectProb = frameObject{j};
    
    % (Appearance + Motion) Probability based object edge detection
%     TargetObject = ObjectProb>0.5;
    
%     [L,numtarget]=bwlabel(TargetObject,8);
%     [m1, n1] = size(ObjectProb);
%     % Delect and clean noisy pixels
%     stats = regionprops(L,'Area');  % Compute the size of each connected region
%     area = cat(1,stats.Area);     
%     index = find(area > 100);       % Finding the index of the region that its size larger than a threshold (e.g.10*10)   
%     % Obtaining the indexed regions which eliminate small noisy regions 
%     TargetObjectBW = ismember(L,index(:));
%     [L,numtarget]=bwlabel(TargetObjectBW,8);
%     % Record cleaned logical labelled Object
%     TargetObject = TargetObjectBW;
    
    TargetObject05  = ObjectProb>0.5;
    TargetObject01  = ObjectProb>0.1;
    TargetObject005 = ObjectProb>0.05;
    TargetObject001 = ObjectProb>0.01;
    
    outputObject05 = sprintf('%s/05Object%05d.bmp',outPath, j);    
    imwrite(TargetObject05, outputObject05);
    outputObject01 = sprintf('%s/01Object%05d.bmp',outPath, j);    
    imwrite(TargetObject01, outputObject01);
    outputObject005 = sprintf('%s/005Object%05d.bmp',outPath, j);    
    imwrite(TargetObject005, outputObject005);
    outputObject001 = sprintf('%s/001Object%05d.bmp',outPath, j);    
    imwrite(TargetObject001, outputObject001);
    
end
