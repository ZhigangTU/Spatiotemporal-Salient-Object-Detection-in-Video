function PixelSeedMask = TMotionConsistentAnalysis(MaskBlock1Pixle2, E_blockmask, PixelMask, uflowdata, vflowdata, normflowdata, data, H, W, E_blockmaskSize)

magflow = normflowdata;
Bmask = E_blockmask;
[m,n] = size( Bmask );
% A_hat = zeros( m, n);
% E_hat = zeros( m, n);
% nBW = 4*sqrt(m/12);
% nBH = 3*sqrt(m/12);
nBH = E_blockmaskSize(1);
nBW = E_blockmaskSize(2);

Seedblock = zeros( m, n);

Pmask = MaskBlock1Pixle2;
[mp,np] = size(Pmask);
nPH = H;
nPW = W;
nSize = ceil(nPW/nBW);

PixelSeedMask = zeros(mp, np);

for j = 1:n-1
    if (max(magflow(:,j)) < 1)
        continue;
    end
        
    Bframe  = reshape(Bmask(:,j),nBH,nBW);
    Pframe  = reshape(Pmask(:,j),nPH,nPW);
    uframe  = reshape(uflowdata(:,j),nPH,nPW);
    vframe  = reshape(vflowdata(:,j),nPH,nPW);
    for jbh = 1:nBH,
    for jbw = 1:nBW,
        if Bframe(jbh,jbw)==1   % >0.6 % =1
%             nhStart=(jbh-1)*nSize + 1;
%             nwStart=(jbw-1)*nSize + 1;
%             
%             OneBlockpixelMask = Pframe(nhStart : nhStart + nSize -1, nwStart : nwStart + nSize -1);
%             OneBlockpixelMask = OneBlockpixelMask>0.6;
%             
%             OneBlocku = uframe(nhStart : nhStart + nSize -1, nwStart : nwStart + nSize -1);
%             OneBlockv = vframe(nhStart : nhStart + nSize -1, nwStart : nwStart + nSize -1);
%             OneBlocku = OneBlocku.*OneBlockpixelMask;
%             OneBlockv = OneBlockv.*OneBlockpixelMask;
%             
%             num = sum(OneBlockpixelMask(:));
%             aveU= sum(OneBlocku(:))/num;
%             aveV= sum(OneBlockv(:))/num;
            
            iter = 0;
            TrackingOver = false;
            curF = j;
            curH = jbh;
            curW = jbw;
            TrackinfoMatrix = [];
            while ~TrackingOver       
                iter = iter + 1;
                bFindnext = false;
                SearchF = curF + 1;
                searchH = curH;
                searchW = curW;
                Bframereset  = reshape(Bmask(:,curF),nBH,nBW);
                % Bframereset(curH,curW) = Bframereset(curH,curW) + 1; 
                Bframereset(curH,curW) = 0;
                Bmask(:,curF)=reshape(Bframereset,nBH*nBW,1);
                
                Bframetemp  = reshape(Bmask(:,SearchF),nBH,nBW);
                Pframecur   = reshape(Pmask(:,curF),nPH,nPW);
                uframecur   = reshape(uflowdata(:,curF),nPH,nPW);
                vframecur   = reshape(vflowdata(:,curF),nPH,nPW);
                [aveU, aveV] = gzComputeAveUV(Pframecur,uframecur,vframecur,curH,curW,nSize); 
                TrackinfoMatrix = [TrackinfoMatrix; curF curH curW aveU aveV];
                xNeibor = [1 1 0 -1 -1 -1 0 1];
                yNeibor = [0 1 1 1 0 -1 -1 -1];
                
                if gzFindMatch( Bframetemp, searchH, searchW )%Bframetemp(searchH,searchW)
                    curF = SearchF;
                    curH = searchH;
                    curW = searchW;
                    bFindnext = 1;
                else
                    nDirection=gzTellDirection(aveU, aveV);
                    nSearch1 = floor(nDirection/10);
                    nSearch2 = mod  (nDirection,10);
%                     % just for debug
%                     if nSearch2==0
%                         [aveU aveV] = gzComputeAveUV(Pframecur, uframecur, vframecur,curH, curW, nSize); 
%                         nDirection=gzTellDirection(aveU, aveV);
%                     end
                    
                    if nSearch1==0
                        searchH1 = curH + 0;
                        searchW1 = curW + 0;
                    else
                        searchH1 = curH + yNeibor(nSearch1);
                        searchW1 = curW + xNeibor(nSearch1);
                    end 
%                     searchH1 = curH + yNeibor(nSearch1);
%                     searchW1 = curW + xNeibor(nSearch1);
                    searchH2 = curH + yNeibor(nSearch2);
                    searchW2 = curW + xNeibor(nSearch2);
                    
                    if gzFindMatch( Bframetemp, searchH1, searchW1 ) %Bframetemp(searchH1,searchW1)
                        curF = SearchF;
                        curH = searchH1;
                        curW = searchW1;
                        bFindnext = 1;
                    else
                        if gzFindMatch( Bframetemp, searchH2, searchW2 ) %Bframetemp(searchH2,searchW2)
                            curF = SearchF;
                            curH = searchH2;
                            curW = searchW2;
                            bFindnext = 1;
                        end
                    end     
                end 
                       
                if ~bFindnext
                    TrackingOver = 1 ; 
                end
                if curF > n-1
                    TrackingOver = 1;
                end 
            end 
            % to analyze the TrackinfoMatrix, 
            % to analyse the motionconsistency
            
            [bConsistent,MaxDiffangle] = gzAnalyseTrackinginfoMatrix(TrackinfoMatrix, n);
            if bConsistent
                % to mark those blocks
                [mMark,nMark] = size( TrackinfoMatrix );
                for jMark = 1:mMark,
                    fMark=TrackinfoMatrix(jMark,1);
                    hMark=TrackinfoMatrix(jMark,2);
                    wMark=TrackinfoMatrix(jMark,3);
                        
                    Seedframe   = reshape( Seedblock(:,fMark), nBH, nBW );
                    Seedframe(hMark, wMark) = 1;
                    Seedblock(:,fMark) = reshape( Seedframe, nBH*nBW,1 ); 
                    
                    imagedata   = reshape( data(:,fMark), nPH, nPW );
                    nhMark = (hMark - 1)*nSize + 1;
                    nwMark = (wMark - 1)*nSize + 1;
                    imagedata(nhMark : nhMark, nwMark : nwMark + nSize -1) = 1;
                    imagedata(nhMark + nSize -1 : nhMark + nSize -1, nwMark : nwMark + nSize -1) = 1;
                    imagedata(nhMark : nhMark + nSize -1, nwMark : nwMark) = 1;
                    imagedata(nhMark : nhMark + nSize -1, nwMark + nSize -1 : nwMark + nSize -1) = 1;
                    
                    data(:,fMark) = reshape( imagedata, nPH*nPW,1 );
                end
                % gaozhi =1;
            end 
            
        end 
            %Blocknorm = norm(BlockTemp,'fro')/nSetBlock;
            %BlockVect = reshape(BlockTemp,nSetBlock*nSetBlock,1);
            %Blocknorm = norm(BlockVect,1)/(nSetBlock*nSetBlock);
            %frameTemp(nhStart : nhStart + nSetBlock-1, nwStart : nwStart + nSetBlock-1)=BlockTemp;
    end
    end
end 

%nSetBlock=8; %the block size will be nSetBlock*nSetBlock
%nWblock=nW/nSetBlock;
%nHblock=nH/nSetBlock;

%%
for j = 1:n,
    Seedframe   = reshape( Seedblock(:,j), nBH, nBW );
    Pixelframe  = reshape( PixelMask(:,j), nPH, nPW );
    
    % treeman case, should clean,
    % CleanSeed = gzCleanSeedblock(Seedframe);
    
    %raincar case do not clean
    CleanSeed   = Seedframe;
    
    cleanBlock  = zeros(nPH, nPW);
    for jh = 1:nBH,
        for jw = 1:nBW,
            nhStart=(jh-1)*nSize+1;
            nwStart=(jw-1)*nSize+1;
            cleanBlock(nhStart:min(nhStart+nSize-1,nPH),nwStart:min(nwStart+nSize-1,nPW))=CleanSeed(jh,jw)*ones(min(nhStart+nSize-1,nPH)-nhStart+1,min(nwStart+nSize-1,nPW)-nwStart+1); %(nSize,nSize)  
        end
    end
    Pixelseed = Pixelframe.*cleanBlock;
%     outputSeed  = sprintf('zzpixelseed%05d.bmp',j);
%     imwrite(Pixelseed, outputSeed);
    
    PixelSeedMask(:,j) = reshape(Pixelseed, nPH*nPW, 1);
    
    % to grow by using the pixel seed
%     TargetPixels = gzSeedgrow(Pixelseed,Pixelframe);
%     outputTarget  = sprintf('ztarget%05d.bmp',j);
%     imwrite(TargetPixels, outputTarget);
    
%     outputMark  = sprintf('zzcleanseed%05d.bmp',j);
%     imwrite(cleanBlock, outputMark);
end

% for j = 1:n,
%   imagedata   = reshape( data(:,j), nPH, nPW );
%   outputMark  = sprintf('wwmotiontrack%05d.bmp',j);
%   imwrite(imagedata, outputMark);
% end

end
