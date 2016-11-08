function CleanSeed = gzCleanSeedblock( Seedframe )

[m n] = size( Seedframe );

CleanSeed=zeros( m, n);

xNeibor=[1 1 0 -1 -1 -1 0 1];
yNeibor=[0 1 1 1 0 -1 -1 -1];
                
for h = 1:m,
    for w = 1:n,
        if Seedframe( h, w )
            nNeibors = 0;
            for k = 1:8,
                hcur = h + yNeibor(k);
                wcur = w + xNeibor(k);
                if hcur<=m&hcur>=1&wcur<=n&wcur>=1
                    if Seedframe( hcur, wcur )
                        nNeibors = nNeibors +1;
                    end 
                end
            end
            if nNeibors==0
                Seedframe( h, w )=0;
            end
        end
    end
end

CleanSeed = Seedframe;
               
end