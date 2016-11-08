% Gao Zhi
% divide and conquer strategy.
% divide original frame into several sub-regions, 
% then, perform rpca for color sequence, to obatin background
% robust batch image alignment example

%for partnum = 0:8
for partnum = 0:11
    % clear
    % clc ; clear all; close all ;
    % addpath
    addpath data ;
    addpath results ;

    %
    methodchoose =5001; %original RPCA code

    %define images' path
    currentPath = cd;
    % input path
    imagePath = fullfile(currentPath,'data') ;

    %userName = 'artscanteen0418ready';
    %userName = '0502rainybeforelunchpeakready';
    %userName = '0506sundayready';
    %userName = 'right1closecollectionready';
    userName = 'engingcam40615';
    
    

    % output path
    if methodchoose==5001
        tmpName  = sprintf('sub.%02dpart',partnum);
        destRoot = fullfile(currentPath,tmpName);
    end
    destDir = fullfile(destRoot,userName) ;
    if ~exist(destDir,'dir')
        mkdir(destRoot,userName) ;
    end
    
    % Get training images
    [fileNames, numImages] = gzget_training_images( imagePath, userName) ;
    
    %% read every images
%     testImage = imread(fileNames{1});
%     if isrgb(testImage)
%         testImage = rgb2gray(testImage);
%     end
%     [h, w] = size(testImage);
    
    %hard code the patch
    if partnum==0
        h=62; w=37;
        hstart = 1; hend=62;
        wstart = 284; wend=320;
    end
    
    if partnum==1
        h=89; w=80;
        hstart = 63; hend=151;
        wstart = 1; wend=80;
    end
    
    if partnum==2
        h=89; w=80;
        hstart = 63; hend=151;
        wstart = 81; wend=160;
    end
    
    if partnum==3
        h=89; w=80;
        hstart = 63; hend=151;
        wstart = 161; wend=240;
    end
    
    if partnum==4
        h=89; w=80;
        hstart = 63; hend=151;
        wstart = 241; wend=320;
    end
    
    if partnum==5
        h=89; w=80;
        hstart = 152; hend=240;
        wstart = 1; wend=80;
    end
    
    if partnum==6
        h=89; w=80;
        hstart = 152; hend=240;
        wstart = 81; wend=160;
    end
    
    if partnum==7
        h=89; w=80;
        hstart = 152; hend=240;
        wstart = 161; wend=240;
    end
    
    if partnum==8
        h=89; w=80;
        hstart = 152; hend=240;
        wstart = 241; wend=320;
    end
    
    if partnum==9
        h=62; w=103;
        hstart = 1; hend=62;
        wstart = 181; wend=283;
    end
    
    if partnum==10
        h=62; w=90;
        hstart = 1; hend=62;
        wstart = 91; wend=180;
    end
    
    if partnum==11
        h=62; w=90;
        hstart = 1; hend=62;
        wstart = 1; wend=90;
    end
    
    dataR = zeros([h*w, numImages], 'double');
    dataG = zeros([h*w, numImages], 'double');
    dataB = zeros([h*w, numImages], 'double');

    for fileIndex = 1:numImages
        currentImage = imread(fileNames{fileIndex});

        %currentImage = currentImage(1:62,284:320,:);
        %currentImage = currentImage(63:151,241:320,:);
        %currentImage = currentImage(152:240,241:320,:);
        currentImage = currentImage(hstart:hend,wstart:wend,:);

        currentImageR = currentImage(:,:,1);
        currentImageG = currentImage(:,:,2);
        currentImageB = currentImage(:,:,3);

        imR = reshape(currentImageR,w*h,1);
        dataR(:, fileIndex) = imR(:,1);

        imG = reshape(currentImageG,w*h,1);
        dataG(:, fileIndex) = imG(:,1);

        imB = reshape(currentImageB,w*h,1);
        dataB(:, fileIndex) = imB(:,1);
    end
    
    dataR=dataR/255.0;
    dataG=dataG/255.0;
    dataB=dataB/255.0;
    
    if methodchoose==5001
    [A_hatR E_hatR iterR] = inexact_alm_rpca(dataR);%original ALM method
    clear iterR;
    
    [A_hatG E_hatG iterG] = inexact_alm_rpca(dataG);%original ALM method
    clear iterG;
    
    [A_hatB E_hatB iterB] = inexact_alm_rpca(dataB);%original ALM method
    %clear E_hatB;
    clear iterB;
    %clear A_hatB  %clear dataB;
    end
    
    dRatio = 1.0;%%0.5;%%
    E_hatR = abs(E_hatR);
    meanR  = mean2(E_hatR);
    stdR   = std2(E_hatR);
    ThreshR=meanR+dRatio*stdR;
    E_hatR=(E_hatR<ThreshR);

    E_hatG = abs(E_hatG);
    meanG  = mean2(E_hatG);
    stdG   = std2(E_hatG);
    ThreshG=meanG+dRatio*stdG;
    E_hatG=(E_hatG<ThreshG);

    E_hatB = abs(E_hatB);
    meanB  = mean2(E_hatB);
    stdB   = std2(E_hatB);
    ThreshB=meanB+dRatio*stdB;
    E_hatB=(E_hatB<ThreshB);

    MaskSafe=E_hatR.*E_hatG;
    MaskSafe=MaskSafe.*E_hatB;

    MaskSafeAnti= 1-MaskSafe;


    dataR =dataR.*MaskSafe+MaskSafeAnti;
    dataG =dataG.*MaskSafe+MaskSafeAnti;
    dataB =dataB.*MaskSafe+MaskSafeAnti;

    WithMask(:,:,1)=dataR* 255;
    WithMask(:,:,2)=dataG* 255;
    WithMask(:,:,3)=dataB* 255;
    MatrixMASKName  = sprintf('WithMask%02d.mat',partnum);
    %save('WithMask.mat','WithMask') ;
    save(MatrixMASKName,'WithMask') ;
    clear WithMask;

    ClearBK(:,:,1)=A_hatR* 255;
    ClearBK(:,:,2)=A_hatG* 255;
    ClearBK(:,:,3)=A_hatB* 255;
    MatrixBKName  = sprintf('Cleanbk%02d.mat',partnum);
    %save('Cleanbk.mat','ClearBK') ;
    save(MatrixBKName,'ClearBK') ;
    clear ClearBK;
    
    %%save the data as images (low rank and sparse matrix)
    for fileIndex = 1:numImages
        imRbackground = A_hatR(:,fileIndex);
        imRbackground = reshape(imRbackground,h,w);
        imRbackground = imRbackground * 255;
        imR = dataR(:,fileIndex);
        imR = reshape(imR,h,w);
        imR = imR * 255;

        imGbackground = A_hatG(:,fileIndex);
        imGbackground = reshape(imGbackground,h,w);
        imGbackground = imGbackground * 255;
        imG = dataG(:,fileIndex);
        imG = reshape(imG,h,w);
        imG = imG * 255;

        imBbackground = A_hatB(:,fileIndex);
        imBbackground = reshape(imBbackground,h,w);
        imBbackground = imBbackground * 255;
        imB = dataB(:,fileIndex);
        imB = reshape(imB,h,w);
        imB = imB * 255;

        BkImage(:,:,1)=imR;
        BkImage(:,:,2)=imG;
        BkImage(:,:,3)=imB;
        LowrankImage(:,:,1)=imRbackground;
        LowrankImage(:,:,2)=imGbackground;
        LowrankImage(:,:,3)=imBbackground;

        outputFileNameSparse  = sprintf('Trainingarea%06d.jpg',fileIndex);
        outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);
        %imwrite(im2, outputFileNamesSparse);
        %imwrite(uint8(im2), outputFileNamesSparse);
        imwrite(uint8(BkImage), outputFileNamesSparse);


        outputFileNameLowrank  = sprintf('Lowrank%06d.jpg',fileIndex);
        outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank);
        %imwrite(im2, outputFileNamesSparse);
        %imwrite(uint8(im2), outputFileNamesSparse);
        imwrite(uint8(LowrankImage), outputFileNamesLowrank);
    end
    
    save('partnum.mat','partnum'); 
    clc ; clear all; close all ;
    load('partnum.mat');
end   

