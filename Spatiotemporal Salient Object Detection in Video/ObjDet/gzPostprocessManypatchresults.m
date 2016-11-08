% GAO ZHI merge the resulting mask matrix and cleanbk matrix 
% from 'patch color RPCA'

%% clear
clc ; clear all; close all ;

%% addpath
addpath data ;
addpath results ;

%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;

%userName = '0502rainybeforelunchpeakready';
%userName = '0506sundayready';
%userName = '0502rainydaypartofallready';
%userName = 'right1closecollectionready';
userName = 'engingcam40615';



methodchoose =3001, %original RPCA code

% output path
if methodchoose==3001
    destRoot  = fullfile(currentPath,'PostMergePatchRPCAresults.lowrankpart') ;
    destRoot1 = fullfile(currentPath,'PostMergePatchRPCAresults.withmaskpart') ;
end   

destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end
destDir1 = fullfile(destRoot1,userName) ;
if ~exist(destDir1,'dir')
    mkdir(destRoot1,userName) ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

%% read every images
testImage = imread(fileNames{1});

if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

for fileIndex = 1:numImages
%for fileIndex = 1:50
    imPrepare= zeros(h,w,3, 'double');
    imPrepare(1:62,1:283,:) = 255;
    outputFileName  = sprintf('gztestlowrank%06d.jpg',fileIndex);
    outputFileNames= fullfile(destDir, outputFileName);     
    imwrite(uint8(imPrepare), outputFileNames);
    
%     imPrepareMask= zeros(h,w,3, 'double');
%     imPrepareMask(1:62,1:283,:) = 255;
    outputFileName  = sprintf('gztestwithmask%06d.jpg',fileIndex);
    outputFileNames= fullfile(destDir1, outputFileName);     
    imwrite(uint8(imPrepare), outputFileNames);
end

%%
for partnum = 0:11
    %clc ; clear all; close all ;
    currentPath = cd;
    %userName = '0506sundayready';
    imagePath = fullfile(currentPath,'PostMergePatchRPCAresults.lowrankpart') ;
    [fileNames, numImages] = gzget_training_images( imagePath, userName) ;
    
    currentPath = cd;
    imagePath1 = fullfile(currentPath,'PostMergePatchRPCAresults.withmaskpart') ;
    [fileNames1, numImages1] = gzget_training_images( imagePath1, userName) ;
    

    FileNameMask  = sprintf('WithMask%02d.mat',partnum);
    %load('WithMask.mat'); %save('WithMask.mat','WithMask') ;
    load(FileNameMask);

    FileNameBK  = sprintf('Cleanbk%02d.mat',partnum);
    % save('Cleanbk.mat','ClearBK') ;
    load(FileNameBK);

    %get two matrix: WithMask, and ClearBK.
    htemp=89;
    wtemp=80;
    if partnum==0
        htemp=62;
        wtemp=37;
    end
    
    if partnum==9
        htemp=62; wtemp=103;
    end
    if partnum==10
        htemp=62;wtemp=90;
    end
    if partnum==11
        htemp=62;wtemp=90;
    end
    
    for fileIndex = 1:numImages
        currentImage = imread(fileNames{fileIndex});
        currentImageMask = imread(fileNames1{fileIndex});

        tempMask = WithMask(:,fileIndex,1);
        tempBK   = ClearBK(:,fileIndex,1);
        tempMask = reshape(tempMask,htemp,wtemp);
        tempBK = reshape(tempBK,htemp,wtemp);
        if partnum==0
            currentImage(1:62,284:320,1) = tempBK;
            currentImageMask(1:62,284:320,1) = tempMask;
        end 
        
        if partnum==1
            currentImage(63:151,1:80,1) = tempBK;
            currentImageMask(63:151,1:80,1) = tempMask;
        end 
        
        if partnum==2
            currentImage(63:151,81:160,1) = tempBK;
            currentImageMask(63:151,81:160,1) = tempMask;
        end 
        
        if partnum==3
            currentImage(63:151,161:240,1) = tempBK;
            currentImageMask(63:151,161:240,1) = tempMask;
        end 
        
        if partnum==4
            currentImage(63:151,241:320,1) = tempBK;
            currentImageMask(63:151,241:320,1) = tempMask;
        end 
        
        if partnum==5
            currentImage(152:240,1:80,1) = tempBK;
            currentImageMask(152:240,1:80,1) = tempMask;
        end 
        
        if partnum==6
            currentImage(152:240,81:160,1) = tempBK;
            currentImageMask(152:240,81:160,1) = tempMask;
        end 
        
        if partnum==7
            currentImage(152:240,161:240,1) = tempBK;
            currentImageMask(152:240,161:240,1) = tempMask;
        end 
        
        if partnum==8
            currentImage(152:240,241:320,1) = tempBK;
            currentImageMask(152:240,241:320,1) = tempMask;
        end 
        
        if partnum==9
            htemp=62; wtemp=103;
            currentImage(1:62,181:283,1) = tempBK;
            currentImageMask(1:62,181:283,1) = tempMask;
        end 
        
        if partnum==10
            htemp=62;wtemp=90;
            currentImage(1:62,91:180,1) = tempBK;
            currentImageMask(1:62,91:180,1) = tempMask;
        end 
        
        if partnum==11
            htemp=62;wtemp=90;
            currentImage(1:62,1:90,1) = tempBK;
            currentImageMask(1:62,1:90,1) = tempMask;
        end 
        
        %
        tempMask = WithMask(:,fileIndex,2);
        tempBK   = ClearBK(:,fileIndex,2);
        tempMask = reshape(tempMask,htemp,wtemp);
        tempBK = reshape(tempBK,htemp,wtemp);
        if partnum==0
            currentImage(1:62,284:320,2) = tempBK;
            currentImageMask(1:62,284:320,2) = tempMask;
        end
        
        if partnum==1
            currentImage(63:151,1:80,2) = tempBK;
            currentImageMask(63:151,1:80,2) = tempMask;
        end
        
        if partnum==2
            currentImage(63:151,81:160,2) = tempBK;
            currentImageMask(63:151,81:160,2) = tempMask;
        end 
        
        if partnum==3
            currentImage(63:151,161:240,2) = tempBK;
            currentImageMask(63:151,161:240,2) = tempMask;
        end 
        
        if partnum==4
            currentImage(63:151,241:320,2) = tempBK;
            currentImageMask(63:151,241:320,2) = tempMask;
        end 
        
        if partnum==5
            currentImage(152:240,1:80,2) = tempBK;
            currentImageMask(152:240,1:80,2) = tempMask;
        end 
        
        if partnum==6
            currentImage(152:240,81:160,2) = tempBK;
            currentImageMask(152:240,81:160,2) = tempMask;
        end 
        
        if partnum==7
            currentImage(152:240,161:240,2) = tempBK;
            currentImageMask(152:240,161:240,2) = tempMask;
        end 
        
        if partnum==8
            currentImage(152:240,241:320,2) = tempBK;
            currentImageMask(152:240,241:320,2) = tempMask;
        end
        
        if partnum==9
            htemp=62; wtemp=103;
            currentImage(1:62,181:283,2) = tempBK;
            currentImageMask(1:62,181:283,2) = tempMask;
        end 
        
        if partnum==10
            htemp=62;wtemp=90;
            currentImage(1:62,91:180,2) = tempBK;
            currentImageMask(1:62,91:180,2) = tempMask;
        end 
        
        if partnum==11
            htemp=62;wtemp=90;
            currentImage(1:62,1:90,2) = tempBK;
            currentImageMask(1:62,1:90,2) = tempMask;
        end
        
        %
        tempMask = WithMask(:,fileIndex,3);
        tempBK   = ClearBK(:,fileIndex,3);
        tempMask = reshape(tempMask,htemp,wtemp);
        tempBK = reshape(tempBK,htemp,wtemp);
        if partnum==0
            currentImage(1:62,284:320,3) = tempBK;
            currentImageMask(1:62,284:320,3) = tempMask;
        end
        
        if partnum==1
            currentImage(63:151,1:80,3) = tempBK;
            currentImageMask(63:151,1:80,3) = tempMask;
        end
        
        if partnum==2
            currentImage(63:151,81:160,3) = tempBK;
            currentImageMask(63:151,81:160,3) = tempMask;
        end 
        
        if partnum==3
            currentImage(63:151,161:240,3) = tempBK;
            currentImageMask(63:151,161:240,3) = tempMask;
        end 
        if partnum==4
            currentImage(63:151,241:320,3) = tempBK;
            currentImageMask(63:151,241:320,3) = tempMask;
        end
        
        if partnum==5
            currentImage(152:240,1:80,3) = tempBK;
            currentImageMask(152:240,1:80,3) = tempMask;
        end 
        
        if partnum==6
            currentImage(152:240,81:160,3) = tempBK;
            currentImageMask(152:240,81:160,3) = tempMask;
        end 
        
        if partnum==7
            currentImage(152:240,161:240,3) = tempBK;
            currentImageMask(152:240,161:240,3) = tempMask;
        end 
        
        if partnum==8
            currentImage(152:240,241:320,3) = tempBK;
            currentImageMask(152:240,241:320,3) = tempMask;
        end 
        
        if partnum==9
            htemp=62; wtemp=103;
            currentImage(1:62,181:283,3) = tempBK;
            currentImageMask(1:62,181:283,3) = tempMask;
        end 
        
        if partnum==10
            htemp=62;wtemp=90;
            currentImage(1:62,91:180,3) = tempBK;
            currentImageMask(1:62,91:180,3) = tempMask;
        end 
        
        if partnum==11
            htemp=62;wtemp=90;
            currentImage(1:62,1:90,3) = tempBK;
            currentImageMask(1:62,1:90,3) = tempMask;
        end

        outputFileNameSparse  = sprintf('gztestlowrank%06d.jpg',fileIndex);
        outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);
        imwrite(uint8(currentImage), outputFileNamesSparse);
        
        outputFileNameSparse  = sprintf('gztestwithmask%06d.jpg',fileIndex);
        outputFileNamesSparse= fullfile(destDir1, outputFileNameSparse);
        imwrite(uint8(currentImageMask), outputFileNamesSparse);
        
        

%         outputFileName  = sprintf('gztestlowrank%06d.jpg',fileIndex);
%         outputFileNames= fullfile(destDir, outputFileName);
%         imwrite(uint8(imPrepare), outputFileNames);
%         outputFileName  = sprintf('gztestwithmask%06d.jpg',fileIndex);
%         outputFileNames= fullfile(destDir1, outputFileName);
%         imwrite(uint8(imPrepare), outputFileNames);

    end
    clear WithMask; 
    clear ClearBK;
end

%     imRbackgroundR = zeros([h, w], 'double');
%     imRbackgroundR(1:62,1:283) = 255;

%     imRbackgroundG = zeros([h, w], 'double');
%     imRbackgroundG(1:62,1:283) = 255;
%     
%     imRbackgroundB = zeros([h, w], 'double');
%     imRbackgroundB(1:62,1:283) = 255;


%imwrite(im2, outputFileNamesSparse);
    %imwrite(uint8(im2), outputFileNamesSparse);
    
    
%         if partnum==2
%             tempMask = WithMask2(:,fileIndex);
%             tempBK   = ClearBK2(:,fileIndex);  
%         end
%         
%         if partnum==3
%             tempMask = WithMask3(:,fileIndex);
%             tempBK   = ClearBK3(:,fileIndex);  
%         end
%         
%         if partnum==4
%             tempMask = WithMask4(:,fileIndex);
%             tempBK   = ClearBK4(:,fileIndex);  
%         end
% 
%         if partnum==5
%             tempMask = WithMask5(:,fileIndex);
%             tempBK   = ClearBK5(:,fileIndex);  
%         end
%         
%         if partnum==6
%             tempMask = WithMask6(:,fileIndex);
%             tempBK   = ClearBK6(:,fileIndex);  
%         end
%         
%         if partnum==7
%             tempMask = WithMask7(:,fileIndex);
%             tempBK   = ClearBK7(:,fileIndex);  
%         end
%         
%         if partnum==8
%             tempMask = WithMask8(:,fileIndex);
%             tempBK   = ClearBK8(:,fileIndex);  
%         end


%     imRbackground = CompletedataMASKR(:,fileIndex);
%     imRbackground = reshape(imRbackground,h,w);
%     imR = CompletedataBKR(:,fileIndex);
%     imR = reshape(imR,h,w);
% 
%     imGbackground = CompletedataMASKG(:,fileIndex);
%     imGbackground = reshape(imGbackground,h,w);
%     imG = CompletedataBKG(:,fileIndex);
%     imG = reshape(imG,h,w);
%   
%     imBbackground = CompletedataMASKB(:,fileIndex);
%     imBbackground = reshape(imBbackground,h,w);
%     imB = CompletedataBKB(:,fileIndex);
%     imB = reshape(imB,h,w);
    
    
%     BkImage(:,:,1)=imRbackgroundR;
%     BkImage(:,:,2)=imRbackgroundG;
%     BkImage(:,:,3)=imRbackgroundB;
    
%     LowrankImage(:,:,1)=imRbackground;
%     LowrankImage(:,:,2)=imGbackground;
%     LowrankImage(:,:,3)=imBbackground;
%     outputFileNameSparse  = sprintf('gztestbk%06d.jpg',fileIndex);
% %     outputFileNameSparse  = sprintf('Lowrank%d.jpg',fileIndex);
%     outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
%     %imwrite(im2, outputFileNamesSparse);
%     %imwrite(uint8(im2), outputFileNamesSparse);
%     imwrite(uint8(BkImage), outputFileNamesSparse);
    
    
%     outputFileNameLowrank  = sprintf('WithMask%d.jpg',fileIndex);
%     outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
%     %imwrite(im2, outputFileNamesSparse);
%     %imwrite(uint8(im2), outputFileNamesSparse);
%     imwrite(uint8(LowrankImage), outputFileNamesLowrank);
%   end