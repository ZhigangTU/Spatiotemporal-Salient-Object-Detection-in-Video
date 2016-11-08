function TargetNoInside = gzDeleteInsideblock(TargetTempSave)

TargetNoInside = [];
[m,n] = size(TargetTempSave);
for nT = 1:m,
     blocktemp = TargetTempSave(nT,:);
     minH = blocktemp(1);
     minW = blocktemp(2);
     maxH = blocktemp(3);
     maxW = blocktemp(4);
     bInside = false;
     for nS = 1:m,
         if nS ~= nT
             blockNext = TargetTempSave(nS,:);
             minHnext = blockNext(1);
             minWnext = blockNext(2);
             maxHnext = blockNext(3);
             maxWnext = blockNext(4);
             if (minH>=minHnext) && (minW>=minWnext) && (maxH<=maxHnext) && (maxW<=maxWnext)
                 bInside = true;
                 break;
             end
         end
     end
     
     if ~bInside 
         TargetNoInside = [TargetNoInside; blocktemp];
     end 
end 

end