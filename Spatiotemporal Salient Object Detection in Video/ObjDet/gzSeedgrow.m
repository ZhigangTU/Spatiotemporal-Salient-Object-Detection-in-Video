function TargetPixels = gzSeedgrow(Pixelseed, Pixelframe)

[m, n] = size(Pixelseed);
TargetPixels = zeros(m, n);

xNeibor=[1 1 0 -1 -1 -1 0 1];
yNeibor=[0 1 1 1 0 -1 -1 -1];

%%
changed = true;

while changed       
    changed = false;
 
    for h = 1:m,
    for w = 1:n,
        if Pixelseed(h,w)>0.6 %==1
            for k = 1:8,
                hcur = h + yNeibor(k);
                wcur = w + xNeibor(k);
                if hcur<=m && hcur>=1 && wcur<=n && wcur>=1
                    if Pixelframe( hcur, wcur )>0.9 && Pixelseed( hcur, wcur )<0.1
                        Pixelseed(hcur, wcur) = 1;
                        changed = true;
                        % k=9;   % no use
                        break;
                    end 
                end
            end
            if changed
                break;
            end
        end
    end
    if changed
        break;
    end
    end
  
    %% stop Criterion    
end
%%

% for h = 1:m,
%     for w = 1:n,
%         if Pixelseed(h,w)>0.6 %==1
%             nAnychange = 0;
%             for k = 1:8,
%                 hcur = h + yNeibor(k);
%                 wcur = w + xNeibor(k);
%                 if hcur<=m&hcur>=1&wcur<=n&wcur>=1
%                     if Pixelframe( hcur, wcur )>0.9 & Pixelseed( hcur, wcur )<0.1
%                         nAnychange = 1;
%                         Pixelseed(hcur, wcur) = 1;
%                         %k=9;%no use
%                         break;
%                     end 
%                 end
%             end
%             if nAnychange ==1
%                 h = 1;
%                 w = 1;
%                 break;
%                 %TargetPixels = gzSeedgrow(Pixelseed, Pixelframe)
%             end
%         end
%     end
% end

TargetPixels = Pixelseed;

%%

end