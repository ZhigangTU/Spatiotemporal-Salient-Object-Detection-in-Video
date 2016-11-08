function [minH minW maxH maxW]  = gzGetCoordofTarget(TargetPixels)

[m n] = size(TargetPixels);
minH = m;
minW = n;
maxH = 0;
maxW = 0;

for h = 1:m,
    for w = 1:n,
        if TargetPixels(h,w)>0.6 
            if h>maxH
                maxH = h;
            end
            if h<minH
                minH = h;
            end
            if w>maxW
                maxW = w;
            end
            if w<minW
                minW = w;
            end
        end
    end
end

nExtend=2;
minH = max(1,minH - nExtend);
minW = max(1,minW - nExtend);
maxH = min(m,maxH + nExtend);
maxW = min(n,maxW + nExtend);
 
end