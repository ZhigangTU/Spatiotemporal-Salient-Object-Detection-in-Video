% GAO ZHI merge the resulting mask matrix and cleanbk matrix from 'patch color RPCA'

% clear
clc ; clear all; close all ;

% addpath
addpath data ;
addpath results ;
%

methodchoose =2001, %original RPCA code
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;

%userName = '0502rainyafterlunchpeakready';
%userName = '0502rainybeforelunchpeakready';
%userName = '0502rainyafterlunchpeakagainready';
%userName = '04240502rainyafterlunchpeakagainready';
%userName = '0502rainybeforelunchpeakready';
userName = '0506sundayready';


% output path
if methodchoose==2001
    destRoot = fullfile(currentPath,'PatchRPCAresultsoriginalpart') ;
end   

destDir = fullfile(destRoot,userName) ;

if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;

%% read every images
testImage = imread(fileNames{1});

if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

%numImages=50;

% CompletedataMASKR = zeros([h*w, numImages], 'int8');
% CompletedataMASKR(1:62,1:283) = 255;
% CompletedataBKR = zeros([h*w, numImages], 'int8');
% CompletedataBKR(1:62,1:283) = 255;
% 
% CompletedataMASKG = zeros([h*w, numImages], 'int8');
% CompletedataMASKG(1:62,1:283) = 255;
% CompletedataBKG = zeros([h*w, numImages], 'int8');
% CompletedataBKG(1:62,1:283) = 255;
% 
% CompletedataMASKB = zeros([h*w, numImages], 'int8');
% CompletedataMASKB(1:62,1:283) = 255;
% CompletedataBKB = zeros([h*w, numImages], 'int8');
% CompletedataBKB(1:62,1:283) = 255;

%CompletedataBKB = zeros([h*w, numImages], 'double');

for partnum = 0:4
  FileNameMask  = sprintf('WithMask%d.mat',partnum);
  %load('WithMask.mat');
  %save('WithMask.mat','WithMask') ;
  load(FileNameMask);
  
  FileNameBK  = sprintf('Cleanbk%d.mat',partnum);
  % save('Cleanbk.mat','ClearBK') ;
  load(FileNameBK);
  
  if partnum==0
    WithMaskUpright = WithMask;
    ClearBKUpright  = ClearBK;    
  end
  if partnum==1
    WithMask1 = WithMask;
    ClearBK1  = ClearBK;  
  end
  if partnum==2
    WithMask2 = WithMask;
    ClearBK2  = ClearBK;
  end
  if partnum==3
    WithMask3 = WithMask;
    ClearBK3  = ClearBK;
  end
  if partnum==4
    WithMask4 = WithMask;
    ClearBK4  = ClearBK;
  end
  
  if partnum==5
    WithMask5 = WithMask;
    ClearBK5  = ClearBK;
  end
  
  if partnum==6
    WithMask6 = WithMask;
    ClearBK6  = ClearBK;
  end
  
  if partnum==7
    WithMask7 = WithMask;
    ClearBK7  = ClearBK;
  end
  
  if partnum==8
    WithMask8 = WithMask;
    ClearBK8  = ClearBK;
  end
  
  clear WithMask;
  clear ClearBK;
end


gztest=0;

% CompletedataMASKR(1:62,284:320) = WithMask(:,:,1);
%     CompletedataBKR(1:62,284:320) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(1:62,284:320) = WithMask(:,:,2);
%     CompletedataBKG(1:62,284:320) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(1:62,284:320) = WithMask(:,:,3);
%     CompletedataBKB(1:62,284:320) = ClearBK(:,:,3);
%     
%     CompletedataMASKR(63:151,1:80) = WithMask(:,:,1);
%     CompletedataBKR  (63:151,1:80) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(63:151,1:80) = WithMask(:,:,2);
%     CompletedataBKG  (63:151,1:80) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(63:151,1:80) = WithMask(:,:,3);
%     CompletedataBKB  (63:151,1:80) = ClearBK(:,:,3);
%     CompletedataMASKR(63:151,81:160) = WithMask(:,:,1);
%     CompletedataBKR  (63:151,81:160,:) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(63:151,81:160) = WithMask(:,:,2);
%     CompletedataBKG  (63:151,81:160) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(63:151,81:160) = WithMask(:,:,3);
%     CompletedataBKB  (63:151,81:160) = ClearBK(:,:,3);
%     CompletedataMASKR(63:151,161:240) = WithMask(:,:,1);
%     CompletedataBKR  (63:151,161:240) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(63:151,161:240) = WithMask(:,:,2);
%     CompletedataBKG  (63:151,161:240) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(63:151,161:240) = WithMask(:,:,3);
%     CompletedataBKB  (63:151,161:240) = ClearBK(:,:,3);
%     CompletedataMASKR(63:151,241:320) = WithMask(:,:,1);
%     CompletedataBKR  (63:151,241:320) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(63:151,241:320) = WithMask(:,:,2);
%     CompletedataBKG  (63:151,241:320) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(63:151,241:320) = WithMask(:,:,3);
%     CompletedataBKB  (63:151,241:320) = ClearBK(:,:,3);
%     CompletedataMASKR(152:240,1:80) = WithMask(:,:,1);
%     CompletedataBKR  (152:240,1:80) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(152:240,1:80) = WithMask(:,:,2);
%     CompletedataBKG  (152:240,1:80) = ClearBK(:,:,2);
%     
%     CompletedataMASK(152:240,1:80) = WithMask(:,:,3);
%     CompletedataBK  (152:240,1:80) = ClearBK(:,:,3);
%     CompletedataMASKR(152:240,81:160) = WithMask(:,:,1);
%     CompletedataBKR  (152:240,81:160) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(152:240,81:160) = WithMask(:,:,2);
%     CompletedataBKG  (152:240,81:160) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(152:240,81:160) = WithMask(:,:,3);
%     CompletedataBKB  (152:240,81:160) = ClearBK(:,:,3);
%     CompletedataMASKR(152:240,161:240) = WithMask(:,:,1);
%     CompletedataBKR  (152:240,161:240) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(152:240,161:240) = WithMask(:,:,2);
%     CompletedataBKG  (152:240,161:240) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(152:240,161:240) = WithMask(:,:,3);
%     CompletedataBKB  (152:240,161:240) = ClearBK(:,:,3);
%     CompletedataMASKR(152:240,241:320) = WithMask(:,:,1);
%     CompletedataBKR  (152:240,241:320) = ClearBK(:,:,1);
%     
%     CompletedataMASKG(152:240,241:320) = WithMask(:,:,2);
%     CompletedataBKG  (152:240,241:320) = ClearBK(:,:,2);
%     
%     CompletedataMASKB(152:240,241:320) = WithMask(:,:,3);
%     CompletedataBKB  (152:240,241:320) = ClearBK(:,:,3);

 for fileIndex = 1:numImages
 %for fileIndex = 1:50
    
    imRbackgroundR = zeros([h, w], 'double');
    imRbackgroundR(1:62,1:283) = 255;
    
    imRbackgroundG = zeros([h, w], 'double');
    imRbackgroundG(1:62,1:283) = 255;
    
    imRbackgroundB = zeros([h, w], 'double');
    imRbackgroundB(1:62,1:283) = 255;
    
    for Fillin = 0:4
        htemp=89;
        wtemp=80;
        if Fillin==0
            htemp=62;
            wtemp=37;
            %
            tempMask = WithMaskUpright(:,fileIndex,1);
            tempBK   = ClearBKUpright(:,fileIndex,1);
            
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundR(1:62,284:320) = tempBK;
            %
            tempMask = WithMaskUpright(:,fileIndex,2);
            tempBK   = ClearBKUpright(:,fileIndex,2);
            
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundG(1:62,284:320) = tempBK;
            %
            tempMask = WithMaskUpright(:,fileIndex,3);
            tempBK   = ClearBKUpright(:,fileIndex,3);
            
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundB(1:62,284:320) = tempBK;
            
        end
        
        if Fillin==1
            tempMask = WithMask1(:,fileIndex,1);
            tempBK   = ClearBK1(:,fileIndex,1);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundR(63:151,1:80) = tempBK;
            
            tempMask = WithMask1(:,fileIndex,2);
            tempBK   = ClearBK1(:,fileIndex,2);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundG(63:151,1:80) = tempBK;
            
            tempMask = WithMask1(:,fileIndex,3);
            tempBK   = ClearBK1(:,fileIndex,3);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundB(63:151,1:80) = tempBK;
        end
        
        if Fillin==2
            tempMask = WithMask2(:,fileIndex,1);
            tempBK   = ClearBK2(:,fileIndex,1);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundR(63:151,81:160) = tempBK;
            
            tempMask = WithMask2(:,fileIndex,2);
            tempBK   = ClearBK2(:,fileIndex,2);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundG(63:151,81:160) = tempBK;
            
            tempMask = WithMask2(:,fileIndex,3);
            tempBK   = ClearBK2(:,fileIndex,3);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundB(63:151,81:160) = tempBK;
        end
        
        if Fillin==3
            tempMask = WithMask3(:,fileIndex,1);
            tempBK   = ClearBK3(:,fileIndex,1);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundR(63:151,161:240) = tempBK;
            
            tempMask = WithMask3(:,fileIndex,2);
            tempBK   = ClearBK3(:,fileIndex,2);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundG(63:151,161:240) = tempBK;
            
            tempMask = WithMask3(:,fileIndex,3);
            tempBK   = ClearBK3(:,fileIndex,3);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundB(63:151,161:240) = tempBK;
        end
        
        if Fillin==4
            tempMask = WithMask4(:,fileIndex,1);
            tempBK   = ClearBK4(:,fileIndex,1);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundR(63:151,241:320) = tempBK;
            
            tempMask = WithMask4(:,fileIndex,2);
            tempBK   = ClearBK4(:,fileIndex,2);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundG(63:151,241:320) = tempBK;
            
            tempMask = WithMask4(:,fileIndex,3);
            tempBK   = ClearBK4(:,fileIndex,3);  
            tempMask = reshape(tempMask,htemp,wtemp);
            tempBK = reshape(tempBK,htemp,wtemp);
            imRbackgroundB(63:151,241:320) = tempBK;
        end
        
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
    end 
    

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
    
    
    BkImage(:,:,1)=imRbackgroundR;
    BkImage(:,:,2)=imRbackgroundG;
    BkImage(:,:,3)=imRbackgroundB;
    
%     LowrankImage(:,:,1)=imRbackground;
%     LowrankImage(:,:,2)=imGbackground;
%     LowrankImage(:,:,3)=imBbackground;
    outputFileNameSparse  = sprintf('gztestbk%06d.jpg',fileIndex);
%     outputFileNameSparse  = sprintf('Lowrank%d.jpg',fileIndex);
    outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
    %imwrite(im2, outputFileNamesSparse);
    %imwrite(uint8(im2), outputFileNamesSparse);
    imwrite(uint8(BkImage), outputFileNamesSparse);
    
    
%     outputFileNameLowrank  = sprintf('WithMask%d.jpg',fileIndex);
%     outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
%     %imwrite(im2, outputFileNamesSparse);
%     %imwrite(uint8(im2), outputFileNamesSparse);
%     imwrite(uint8(LowrankImage), outputFileNamesLowrank);
  end


%%save the data as images (low rank and sparse matrix)
% for fileIndex = 1:numImages
%     
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
%     
%     
%     BkImage(:,:,1)=imR;
%     BkImage(:,:,2)=imG;
%     BkImage(:,:,3)=imB;
%     
%     LowrankImage(:,:,1)=imRbackground;
%     LowrankImage(:,:,2)=imGbackground;
%     LowrankImage(:,:,3)=imBbackground;
%     
%     outputFileNameSparse  = sprintf('Lowrank%d.jpg',fileIndex);
%     outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
%     %imwrite(im2, outputFileNamesSparse);
%     %imwrite(uint8(im2), outputFileNamesSparse);
%     imwrite(uint8(BkImage), outputFileNamesSparse);
%     
%     
%     outputFileNameLowrank  = sprintf('WithMask%d.jpg',fileIndex);
%     outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
%     %imwrite(im2, outputFileNamesSparse);
%     %imwrite(uint8(im2), outputFileNamesSparse);
%     imwrite(uint8(LowrankImage), outputFileNamesLowrank);
% end


% save('WithMask.mat','WithMask') ;
% ClearBK(:,:,1)=A_hatR* 255;
% ClearBK(:,:,2)=A_hatG* 255;
% ClearBK(:,:,3)=A_hatB* 255;
% save('Cleanbk.mat','ClearBK') ;
% clear ClearBK;

%hard code the patch
% right-up patch
% h=62;
% w=37;
% h=89;
% w=80;

% dataR = zeros([h*w, numImages], 'double');
% dataG = zeros([h*w, numImages], 'double');
% dataB = zeros([h*w, numImages], 'double');