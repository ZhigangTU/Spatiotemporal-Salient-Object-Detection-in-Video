%function gzOUT = CameraTraining (foldername, bIsWorkingday, bSittingarea )


%clean the environment
%save('inputdata.mat','foldername','bIsWorkingday','bSittingarea') ;       
clc ;
clear all;
close all ;
%load('inputdata.mat');

% addpath
addpath training_data ;
addpath training_maskdata ;
% addpath training_results;

currentPath = cd;

% input path
imagePath = fullfile(currentPath,'training_data') ;

imagePathMask = fullfile(currentPath,'training_maskdata') ;

%userName = 'Frontier_Cam01testaccumulate';
% userName = 'Frontier_Cam01date20130311to0315';
% userName = 'Frontier_Cam01testaccumulate38to73';
% userName = 'Frontier_Cam01testaccumulate1to100';
userName = 'Frontier_Cam01testaccumulate1to200';


destRoot = fullfile(currentPath,'training_accumulate_results') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

[fileNamesMask, numImagesMask] = gzget_training_images( imagePathMask, userName) ;

%% read every images
testImage = imread(fileNames{1});
if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

if h ~= 240 | w~=320
    disp('the dimension of image is not correct! Exit from the training!!') ;
    return;
end

%%
MaskMatrix = zeros(h*w, numImagesMask); 
dataR     = zeros([h*w, numImagesMask], 'double');
%dataRdiff = zeros([h*w, numImagesMask], 'double');
dataG     = zeros([h*w, numImagesMask], 'double');
%dataGdiff = zeros([h*w, numImagesMask], 'double');
dataB     = zeros([h*w, numImagesMask], 'double');
%dataBdiff = zeros([h*w, numImagesMask], 'double');

for fileIndex = 1:numImagesMask
    OriImage = imread(fileNames{fileIndex});
    
    OriImageR = OriImage(:,:,1);
    OriImageR = reshape(OriImageR,w*h,1);
    dataR(:, fileIndex) = OriImageR(:,1);
    
    OriImageG = OriImage(:,:,2);
    OriImageG = reshape(OriImageG,w*h,1);
    dataG(:, fileIndex) = OriImageG(:,1);
    
    OriImageB = OriImage(:,:,3);
    OriImageB = reshape(OriImageB,w*h,1);
    dataB(:, fileIndex) = OriImageB(:,1);
    
    MaskImage = imread(fileNamesMask{fileIndex});
    MaskImage = reshape(MaskImage,w*h,1);
    MaskMatrix(:, fileIndex) = MaskImage(:,1);
end

OutlierCheck=[];
for pixel = 1:h*w
    for fileIndex = 1:numImagesMask-1
        MaskValue = MaskMatrix(pixel,fileIndex);
        nCounttmp=1;
        while MaskValue==1&&fileIndex<numImagesMask
            if nCounttmp==1
                nStart=fileIndex;
                Rstart=dataR(pixel,fileIndex);
                Gstart=dataG(pixel,fileIndex);
                Bstart=dataB(pixel,fileIndex);
            end
            nCounttmp   =nCounttmp+1;
            fileIndex    = fileIndex+1;
            MaskValueNext= MaskMatrix(pixel,fileIndex);
            Rnext=dataR(pixel,fileIndex);
            Gnext=dataG(pixel,fileIndex);
            Bnext=dataB(pixel,fileIndex);
            
            Rdiffperentage=abs(Rstart-Rnext)/Rstart;
            Gdiffperentage=abs(Gstart-Gnext)/Gstart;
            Bdiffperentage=abs(Bstart-Bnext)/Bstart;
            Diffvect=[Rdiffperentage Gdiffperentage Bdiffperentage];
            maxDiffpercentage=max(Diffvect);
            
            %AppearanceDiff = abs(Rstart-Rnext) + abs(Gstart-Gnext) + abs(Bstart-Bnext);
            
            if maxDiffpercentage>0.15 || fileIndex>=numImagesMask ||MaskValueNext<1
                nEnd=fileIndex;
                if nEnd-nStart>2
                    Vectsave=[pixel;nStart;nEnd;nEnd-nStart];
                    OutlierCheck=[OutlierCheck Vectsave];
                end
                
                break;
            end 
            MaskValue=MaskValueNext;
        end
        
        
    end
    
end

AccumulateOutlierTmp = strcat(userName, 'Accumulate.mat'); 
save(AccumulateOutlierTmp,'OutlierCheck') ;
%%
DurationRow = OutlierCheck(4,:);
maxDuration=max(DurationRow);
minDuration=min(DurationRow);

Hist     = zeros([1, maxDuration-minDuration+1], 'double');
for Duration = minDuration:maxDuration
    DurIndicator=DurationRow==Duration;
    num=sum(DurIndicator);
    Hist(Duration-minDuration+1)=num;
end

xaxis=minDuration:maxDuration;
yaxis=Hist;

plot(xaxis, yaxis, '-ro');
%%
for nTh = minDuration:maxDuration
    %     nThreshlast1 = 10;
    nThreshlast1 = nTh;
    DurationSelect = find(DurationRow>nThreshlast1);
    [nRtmp, nNumtmp] = size(DurationSelect);
    OutlierSelect  =OutlierCheck(:,DurationSelect);

    posShow = zeros([240, 320], 'double');

    for num = 1:nNumtmp
        IndexPos=OutlierSelect(1,num);

        nCol = floor((IndexPos-1)/240)+1;
        nRow = rem  (IndexPos-1, 240)+1;
        posShow(nRow, nCol)=1;
        
    end
    imshow(posShow);

    outputFileName  = sprintf('poswithThresh%06d.bmp',nThreshlast1);
    % outputFileNameSparse  = sprintf('background%d.jpg',fileIndex);
    % outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);
    imwrite(posShow, outputFileName);
end



gaozhi=100;
% clear all;
%[i,j]=find(a>0.5);
%%


%             ARNames = fullfile(destDir, 'ARmatrix.mat');
%             ARNames = fullfile(destDir,ARNamestmp);
%             save(ARNames,'AR_hat') ;
% for fileIndex = 2:numImagesMask
%     PreCol=dataR(:, fileIndex-1);
%     CurCol=dataR(:, fileIndex);
%     Coldiff=CurCol-PreCol;
%     dataRdiff(:, fileIndex) = Coldiff;
%     
%     PreCol=dataG(:, fileIndex-1);
%     CurCol=dataG(:, fileIndex);
%     Coldiff=CurCol-PreCol;
%     dataGdiff(:, fileIndex) = Coldiff;
%     
%     PreCol=dataB(:, fileIndex-1);
%     CurCol=dataB(:, fileIndex);
%     Coldiff=CurCol-PreCol;
%     dataBdiff(:, fileIndex) = Coldiff;
% end

% nThresh   = 10;
% dataRdiff = abs(dataRdiff);
% dataGdiff = abs(dataGdiff);
% dataBdiff = abs(dataBdiff);
% Diffmatrix= dataRdiff + dataGdiff + dataBdiff;
% DiffWithin= Diffmatrix < nThresh;

% % nSampleinterval = 10;
% % if bIsWorkingday|bSittingarea
% %     nSampleinterval=15;
% % end
% 
% nSampleinterval = 3;
% if bIsWorkingday|bSittingarea
%     nSampleinterval=3;
% end
% 
% %nSampleinterval=max(nSampleinterval, ceil(numImages/900)),
% nSampleinterval=3,
% nNumberused = floor(numImages/nSampleinterval),
% NameList=[];
% for nCurrentRound = 1:nRound
%     for nColor = 1:3
%         %data = zeros([h*w, numImages], 'double');
%         %data = zeros([h*w, nNumberused], 'double');
%         data = zeros([h*w, 600], 'double');
%         nCount = 1;
%         %for fileIndex = 1:nSampleinterval:numImages
%         for fileIndexCount = 1:600
%             fileIndex = (fileIndexCount-1)*nRound+nCurrentRound;
%             if nColor==1
%                 NameList=[NameList;fileNames{fileIndex}];
%             end
%             currentImage = imread(fileNames{fileIndex});
%             currentImage = currentImage(:,:,nColor);
%             im = reshape(currentImage,w*h,1);
%             %data(:, fileIndex) = im(:,1);
%             data(:, nCount) = im(:,1);
%             nCount=nCount+1;
%         end
%         data=data/255.0;
%         
%         
%         if nColor==1
%             [AR_hat ER_hat iterR] = inexact_alm_rpca(data);
%             
%             ARNamestmp = sprintf('ARmatrix%02d.mat',nCurrentRound); 
% %             ARNames = fullfile(destDir, 'ARmatrix.mat');
%             ARNames = fullfile(destDir,ARNamestmp);
%             save(ARNames,'AR_hat') ;
%             
%             ERNamestmp = sprintf('ERmatrix%02d.mat',nCurrentRound); 
% %             ERNames = fullfile(destDir, 'ERmatrix.mat');
%             ERNames = fullfile(destDir, ERNamestmp);
%             save(ERNames,'ER_hat') ;
%             clear AR_hat; clear ER_hat; clear iterR;
%         end
%         
%         if nColor==2
%             [AG_hat EG_hat iterG] = inexact_alm_rpca(data);
%             
%             AGNamestmp = sprintf('AGmatrix%02d.mat',nCurrentRound);
%             AGNames = fullfile(destDir, AGNamestmp);
% %             AGNames = fullfile(destDir, 'AGmatrix.mat');
%             save(AGNames,'AG_hat') ;
%             
%             EGNamestmp = sprintf('EGmatrix%02d.mat',nCurrentRound);
%             EGNames = fullfile(destDir, EGNamestmp);
% %             EGNames = fullfile(destDir, 'EGmatrix.mat');
%             save(EGNames,'EG_hat') ;
%             clear AG_hat; clear EG_hat; clear iterG;
%         end
%         
%         if nColor==3
%             [AB_hat EB_hat iterB] = inexact_alm_rpca(data);
%             
%             ABNamestmp = sprintf('ABmatrix%02d.mat',nCurrentRound);
%             ABNames = fullfile(destDir, ABNamestmp);
% %             ABNames = fullfile(destDir, 'ABmatrix.mat');
%             save(ABNames,'AB_hat') ;
%             
%             EBNamestmp = sprintf('EBmatrix%02d.mat',nCurrentRound);
%             EBNames = fullfile(destDir, EBNamestmp);
% %             EBNames = fullfile(destDir, 'EBmatrix.mat');
%             save(EBNames,'EB_hat') ;
%             clear AB_hat; clear EB_hat; clear iterB;
%         end
%         clear data;
%         
%     end
% end
% 
% for nCurrentRound = 1:nRound
%     %EBfilename = fullfile(dataPath, 'EBmatrix') ;
%     EBNamestmp = sprintf('EBmatrix%02d.mat',nCurrentRound);
%     EBNames = fullfile(destDir, EBNamestmp);
%     %load(EBfilename);
%     load(EBNames);
%     eb = abs(EB_hat);
%     clear EB_hat;
%     
%     EGNamestmp = sprintf('EGmatrix%02d.mat',nCurrentRound);
%     EGNames = fullfile(destDir, EGNamestmp);
%     %EGfilename = fullfile(dataPath, 'EGmatrix') ;
%     %load(EGfilename);
%     load(EGNames);
%     eg = abs(EG_hat);
%     clear EG_hat;
%     
%     ERNamestmp = sprintf('ERmatrix%02d.mat',nCurrentRound);
%     ERNames = fullfile(destDir, ERNamestmp);
%     %     ERfilename = fullfile(dataPath, 'ERmatrix') ;
%     %     load(ERfilename);
%     load(ERNames);
%     er = abs(ER_hat);
%     clear ER_hat;
%     
%     eMat = eb+eg+er;
%     ebthreshold = mean(mean(eMat))*1.2;
%     eoutliermatrix = eMat>ebthreshold;
%     clear eb; clear eg; clear er;
%     
%     MaskNamestmp = sprintf('Maskmatrix%02d.mat',nCurrentRound);
%     MaskEBNames = fullfile(destDir, MaskNamestmp);
%     save(MaskEBNames,'eoutliermatrix') ;
%     
%     %     for fileIndexCount = 1:600
%     for fileIndexCount = 1:nSeqNum
%         maskCurrent = eoutliermatrix(:,fileIndexCount);
%         maskCurrent = reshape(maskCurrent,240,320);
%         
%         fileIndex = (fileIndexCount-1)*nRound+nCurrentRound;
%         
%         Maskfiletmp = sprintf('Mask%06d.bmp',fileIndex);
%         Maskfile    = fullfile(destDir, Maskfiletmp);
%         imwrite(maskCurrent, Maskfile);
%     end
%     
%     clear eoutliermatrix;
%     
% end
% load(ARNames);
% load(ERNames);
% AR_hat=abs(AR_hat);
% ER_hat=abs(ER_hat);
% ERsum = sum(ER_hat);
% 
% load(AGNames);
% load(EGNames);
% AG_hat=abs(AG_hat);
% EG_hat=abs(EG_hat);
% 
% load(ABNames);
% load(EBNames);
% AB_hat=abs(AB_hat);
% EB_hat=abs(EB_hat);
% 
% [hHuge, wHuge] = size(AB_hat);
% 
% %%
% %for further use
% Eenergy = ER_hat + EG_hat + EB_hat;
% fEaverage= mean(Eenergy(:));
% Etarget = Eenergy>fEaverage;
% 
% Erowsum = sum(Etarget,2);
% fTh = nNumberused * 0.08;
% Erowlargeoutlier=Erowsum>fTh;
% Erowlargeoutliersum=sum(Erowlargeoutlier),
% 
% 
% txtname     =  strcat(userName,'number.txt');
% fulltxtname = fullfile(destDir, txtname);
% fid = fopen(fulltxtname','wt'); 
% fprintf(fid,'%d\n',Erowlargeoutliersum); 
% fclose(fid);
% 
% %%
% 
% ERsum = sum(ER_hat);
% EGsum = sum(EG_hat);
% EBsum = sum(EB_hat);
% Esum  = ERsum + EGsum + EBsum;
% 
% Esort = sort(Esum);
% 
% fRatio=0.85;
% if bIsWorkingday|bSittingarea
%     fRatio=0.75;
% end
% 
% hThresh=Esort(floor(wHuge*fRatio));
% 
% %if Ekeep(i)=1, keep for Color RPCA
% Ekeep=Esum<hThresh;
% Keepnum=sum(Ekeep), 
% 
% %%second pass RPCA
% clear AR_hat; clear ER_hat;
% clear AG_hat; clear EG_hat;
% clear AB_hat; clear EB_hat;
% 
% for nColor = 1:3
%     
%     data = zeros([h*w, Keepnum], 'double');
%     nCount = 1;
%     for fileIndex = 1:wHuge
%         if Ekeep(fileIndex)==1
%             currentImage = imread(NameList(fileIndex,:));
%             currentImage = currentImage(:,:,nColor);
%             im = reshape(currentImage,w*h,1);
%             %data(:, fileIndex) = im(:,1);
%             data(:, nCount) = im(:,1);
%             nCount=nCount+1;
%         end   
%     end
%     data=data/255.0;
%     
%     if nColor==1
%         [AR_hat ER_hat iterR] = inexact_alm_rpca(data);
%         
%         ARNames = fullfile(destDir, 'ARmatrix.mat');
%         save(ARNames,'AR_hat') ;
%         ERNames = fullfile(destDir, 'ERmatrix.mat');
%         save(ERNames,'ER_hat') ;
%         clear AR_hat; clear ER_hat; clear iterR;
%     end
%     
%     if nColor==2
%         [AG_hat EG_hat iterG] = inexact_alm_rpca(data);
%         
%         AGNames = fullfile(destDir, 'AGmatrix.mat');
%         save(AGNames,'AG_hat') ;
%         EGNames = fullfile(destDir, 'EGmatrix.mat');
%         save(EGNames,'EG_hat') ;
%         clear AG_hat; clear EG_hat; clear iterG;
%     end
%     
%     if nColor==3
%         [AB_hat EB_hat iterB] = inexact_alm_rpca(data);
%         
%         ABNames = fullfile(destDir, 'ABmatrix.mat');
%         save(ABNames,'AB_hat') ;
%         EBNames = fullfile(destDir, 'EBmatrix.mat');
%         save(EBNames,'EB_hat') ;
%         clear AB_hat; clear EB_hat; clear iterB;
%     end
%     clear data;
%     
% end
% 
% %%%%
% 
% load(ARNames);
% load(ERNames);
% AR_hat=abs(AR_hat);
% ER_hat=abs(ER_hat);
% ERsum = sum(ER_hat);
% 
% load(AGNames);
% load(EGNames);
% AG_hat=abs(AG_hat);
% EG_hat=abs(EG_hat);
% 
% load(ABNames);
% load(EBNames);
% AG_hat=abs(AG_hat);
% EG_hat=abs(EG_hat);
% 
% [hHuge, wHuge] = size(AB_hat);
% 
% 
% for fileIndex = 1:wHuge
% %     imR = ER_hat(:,fileIndex);
% %     imR = reshape(imR,h, w);
% %     
% %     imG = EG_hat(:,fileIndex);
% %     imG = reshape(imG,h, w);
% %     
% %     imB = EB_hat(:,fileIndex);
% %     imB = reshape(imB,h, w);
% %     
% %     imageE(:,:,1)=imR;
% %     imageE(:,:,2)=imG;
% %     imageE(:,:,3)=imB;
% %     
% %     outputFileNameLowrank  = sprintf('Lowrank%d.jpg',fileIndex);
% %     outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
% %     imwrite(imageE, outputFileNamesLowrank);
%     
%     imR = AR_hat(:,fileIndex);
%     imR = reshape(imR,h, w);
%     
%     imG = AG_hat(:,fileIndex);
%     imG = reshape(imG,h, w);
%     
%     imB = AB_hat(:,fileIndex);
%     imB = reshape(imB,h, w);
%     
%     imageBK(:,:,1)=imR;
%     imageBK(:,:,2)=imG;
%     imageBK(:,:,3)=imB;
%     
% %     im2 = A_hat(:,fileIndex);
% %     im2 = reshape(im2,h,w);
%     imageBK = imageBK * 255;
%     
%     outputtemp = sprintf('bk%06d.jpg',fileIndex);
%     
%     outputFileNameSparse =  strcat(userName,outputtemp);
%     outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);
%     %imwrite(im2, outputFileNamesSparse);
%     imwrite(uint8(imageBK), outputFileNamesSparse);
% end
% 
% disp('you successfully finish camera training with data in this folder!!') ;
% %end
