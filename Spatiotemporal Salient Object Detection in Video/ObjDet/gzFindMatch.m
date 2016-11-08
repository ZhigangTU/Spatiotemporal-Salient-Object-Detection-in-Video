function bfind = gzFindMatch( DataMatrix, posH, posW )
%bFindMatch(Bframetemp,searchH,searchW)
[m n] = size( DataMatrix );
bfind = 0;
if posH>m | posH<1 | posW>n | posW<1
    bfind = 0;
else
    if DataMatrix(posH,posW)>=1
        bfind=1;
    end
end

end