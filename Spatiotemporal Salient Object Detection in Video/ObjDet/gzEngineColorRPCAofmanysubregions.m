% Gao Zhi
% For engine canteen
% divide and conquer strategy.
% divide original frame into 4*3 sub-regions, 
% then, perform rpca for color sequence, to obatin background
% robust batch image alignment example

% Rmeanall = zeros([1, numImages*12], 'double');
% Gmeanall = zeros([1, numImages*12], 'double');
% Bmeanall = zeros([1, numImages*12], 'double');

Rmeanall = [];
Gmeanall = [];
Bmeanall = [];
nRegionsAll  = 0;

for partnum = 1:12
    % clear% clc ; clear all; close all ;
    addpath data ;
    methodchoose =5001; %original RPCA code

    %define images' path
    currentPath = cd;
    % input path
    imagePath = fullfile(currentPath,'data') ;

    %userName  = 'test8Junepart';
    %userName  = '13Junepart';
    %userName  = 'cam2on4august';
    userName  = 'cam01on6789august';
    
    
    destRoot = fullfile(currentPath,userName);
    subName  = sprintf('sub.%02dpart',partnum);
    
    destDir = fullfile(destRoot,subName) ;
    if ~exist(destDir,'dir')
        mkdir(destRoot,subName) ;
    end
    
    % Get training images
    [fileNames, numImages] = gzget_training_images( imagePath, userName) ;
    nRegionsAll = numImages*12;
    
    %hard code the patch
    h=80; 
    w=80;
    
    shang = floor((partnum-1)/4);
    yushu = rem(partnum-1,  4);
    
    hstart = shang*h+1; 
    hend   = (shang+1)*h;
    wstart = yushu*w+ 1; 
    wend   = (yushu+1)*w;
    
%     if partnum==1
%         hstart = 1; hend=80;
%         wstart = 1; wend=80;
%     end
%     if partnum==2
%         hstart = 1; hend=80;
%         wstart = 81; wend=160;
%     end
%     if partnum==3
%         hstart = 1; hend=80;
%         wstart = 161; wend=240;
%     end
%     if partnum==4
%         hstart = 1; hend=80;
%         wstart = 241; wend=320;
%     end
%     
%     if partnum==5
%         hstart = 81; hend=160;
%         wstart = 1; wend=80;
%     end
%     if partnum==6
%         hstart = 81; hend=160;
%         wstart = 81; wend=160;
%     end
%     if partnum==7
%         hstart = 81; hend=160;
%         wstart = 161; wend=240;
%     end
%     if partnum==8
%         hstart = 81; hend=160;
%         wstart = 241; wend=320;
%     end
%     
%     if partnum==9
%         hstart = 161; hend=240;
%         wstart = 1; wend=80;
%     end
%     if partnum==10
%         hstart = 161; hend=240;
%         wstart = 81; wend=160;
%     end
%     if partnum==11
%         hstart = 161; hend=240;
%         wstart = 161; wend=240;
%     end
%     if partnum==12
%         hstart = 161; hend=240;
%         wstart = 241; wend=320;
%     end
    dataR = zeros([h*w, numImages], 'double');
    dataG = zeros([h*w, numImages], 'double');
    dataB = zeros([h*w, numImages], 'double');

    for fileIndex = 1:numImages
        currentImage = imread(fileNames{fileIndex});
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
        clear iterR;clear dataR;clear A_hatR;
        %
        [A_hatG E_hatG iterG] = inexact_alm_rpca(dataG);%original ALM method
        clear iterG;clear dataG;clear A_hatG;
        %
        [A_hatB E_hatB iterB] = inexact_alm_rpca(dataB);%original ALM method
        clear iterB;clear dataB;clear A_hatB;
    end
    
    MatrixERname  = sprintf('MatrixER%02d.mat',partnum);
    outputMatrixERname = fullfile(destDir, MatrixERname); 
    save(outputMatrixERname,'E_hatR') ;
    RMean = mean(E_hatR,1);
    Rmeanall = [Rmeanall, RMean];
    
    MatrixEGname  = sprintf('MatrixEG%02d.mat',partnum);
    outputMatrixEGname = fullfile(destDir, MatrixEGname); 
    save(outputMatrixEGname,'E_hatG') ;
    GMean = mean(E_hatG,1);
    Gmeanall = [Gmeanall, GMean];

    MatrixEBname  = sprintf('MatrixEB%02d.mat',partnum);
    outputMatrixEBname = fullfile(destDir, MatrixEBname); 
    save(outputMatrixEBname,'E_hatB') ;
    BMean = mean(E_hatB,1);
    Bmeanall = [Bmeanall, BMean];
    clear E_hatR;
    clear E_hatG;
    clear E_hatB;   
end  

yaxis=[Rmeanall, Gmeanall, Bmeanall];
xaxis=[1:nRegionsAll*3];
%plot(xaxis,yaxis);

yvect = yaxis';
meanvalue = mean(yvect);
stdvalue  = std (yvect);
HighThresh= meanvalue + 0.35*stdvalue;
LowThresh = meanvalue - 0.35*stdvalue;

Check1 = yaxis<HighThresh;
Check2 = yaxis> LowThresh;
Check  = Check1 & Check2;

CheckR =Check(1 : nRegionsAll);
CheckG =Check(nRegionsAll+1 : 2*nRegionsAll);
CheckB =Check(nRegionsAll*2 + 1 : 3*nRegionsAll);
CheckRGB =CheckR & CheckG & CheckB ;

CheckFrame = [];
for partnum = 1:12
    CheckFrame=[CheckFrame;CheckRGB((partnum-1)*numImages +1: partnum*numImages)];
end
Colsum = sum(CheckFrame,1);
FrameOK= Colsum>10;
%FrameOK= Colsum>10;

destRoot = fullfile(currentPath,userName);
destDir = fullfile(destRoot,'afterfiltering') ;
if ~exist(destDir,'dir')
    mkdir(destRoot,'afterfiltering') ;
end

for fileIndex = 1:numImages
    currentImage = imread(fileNames{fileIndex});
    if FrameOK(fileIndex)==1
        outputFileName  = sprintf('afterfiltering%06d.jpg',fileIndex);
        outputFileNames = fullfile(destDir, outputFileName);
        imwrite(currentImage, outputFileNames);
    end
end

%Nonzero=sum(Check');
%numRegions=nRegionsAll*3;
%plot(xaxis,Check);
%Higharray = ones([1, nRegionsAll*3], 'double')*HighThresh;
%Lowarray  = ones([1, nRegionsAll*3], 'double')*LowThresh;
%plot(xaxis,yaxis,'b-',xaxis,Higharray,'r-', xaxis,Lowarray,'r-');
% xaxis=[1:nRegionsAll];
% plot(xaxis,Rmeanall,'r-',xaxis,Gmeanall,'g-', xaxis,Bmeanall,'b-');
gaozhi =100;

%% second pass Color RPCA:
for partnum = 1:12
    % Get training images
    [fileNames2, numImages2] = gzget_training_images( destRoot, 'afterfiltering') ;
    destRoot2 = fullfile(destRoot,'afterfiltering') ;
    
    h=80; 
    w=80;
    if partnum==1
        hstart = 1; hend=80;
        wstart = 1; wend=80;
    end
    if partnum==2
        hstart = 1; hend=80;
        wstart = 81; wend=160;
    end
    if partnum==3
        hstart = 1; hend=80;
        wstart = 161; wend=240;
    end
    if partnum==4
        hstart = 1; hend=80;
        wstart = 241; wend=320;
    end
    
    if partnum==5
        hstart = 81; hend=160;
        wstart = 1; wend=80;
    end
    if partnum==6
        hstart = 81; hend=160;
        wstart = 81; wend=160;
    end
    if partnum==7
        hstart = 81; hend=160;
        wstart = 161; wend=240;
    end
    if partnum==8
        hstart = 81; hend=160;
        wstart = 241; wend=320;
    end
    
    if partnum==9
        hstart = 161; hend=240;
        wstart = 1; wend=80;
    end
    if partnum==10
        hstart = 161; hend=240;
        wstart = 81; wend=160;
    end
    if partnum==11
        hstart = 161; hend=240;
        wstart = 161; wend=240;
    end
    if partnum==12
        hstart = 161; hend=240;
        wstart = 241; wend=320;
    end
    
    dataR2 = zeros([h*w, numImages2], 'double');
    dataG2 = zeros([h*w, numImages2], 'double');
    dataB2 = zeros([h*w, numImages2], 'double');

    for fileIndex = 1:numImages2
        currentImage = imread(fileNames2{fileIndex});
        %currentImage = currentImage(1:62,284:320,:);
        %currentImage = currentImage(63:151,241:320,:);
        %currentImage = currentImage(152:240,241:320,:);
        currentImage = currentImage(hstart:hend,wstart:wend,:);

        currentImageR = currentImage(:,:,1);
        currentImageG = currentImage(:,:,2);
        currentImageB = currentImage(:,:,3);

        imR = reshape(currentImageR,w*h,1);
        dataR2(:, fileIndex) = imR(:,1);
        imG = reshape(currentImageG,w*h,1);
        dataG2(:, fileIndex) = imG(:,1);
        imB = reshape(currentImageB,w*h,1);
        dataB2(:, fileIndex) = imB(:,1);
    end
    dataR=dataR2/255.0;
    dataG=dataG2/255.0;
    dataB=dataB2/255.0;
    
    if methodchoose==5001
        [A_hatR E_hatR iterR] = inexact_alm_rpca(dataR);%original ALM method
        clear iterR;clear dataR;

        [A_hatG E_hatG iterG] = inexact_alm_rpca(dataG);%original ALM method
        clear iterG;clear dataG;

        [A_hatB E_hatB iterB] = inexact_alm_rpca(dataB);%original ALM method
        clear iterB;clear dataB;
    end
  
    MatrixERname  = sprintf('MatrixER%02d.mat',partnum);
    outputMatrixERname = fullfile(destRoot2, MatrixERname); 
    save(outputMatrixERname,'E_hatR') ;
    MatrixARname  = sprintf('MatrixAR%02d.mat',partnum);
    outputMatrixARname = fullfile(destRoot2, MatrixARname); 
    save(outputMatrixARname,'A_hatR') ;
   
    MatrixEGname  = sprintf('MatrixEG%02d.mat',partnum);
    outputMatrixEGname = fullfile(destRoot2, MatrixEGname); 
    save(outputMatrixEGname,'E_hatG') ;
    MatrixAGname  = sprintf('MatrixAG%02d.mat',partnum);
    outputMatrixAGname = fullfile(destRoot2, MatrixAGname); 
    save(outputMatrixAGname,'A_hatG') ;
    
    MatrixEBname  = sprintf('MatrixEB%02d.mat',partnum);
    outputMatrixEBname = fullfile(destRoot2, MatrixEBname); 
    save(outputMatrixEBname,'E_hatB') ;
    MatrixABname  = sprintf('MatrixAB%02d.mat',partnum);
    outputMatrixABname = fullfile(destRoot2, MatrixABname); 
    save(outputMatrixABname,'A_hatB') ;
    
    clear A_hatR;
    clear E_hatR;
    clear A_hatG;
    clear E_hatG;
    clear A_hatB;
    clear E_hatB;
end  

%% after 2-pass COLOR RPCA, 
%  post processing 
destDirlowrank = fullfile(destRoot,'PostMergePatchRPCAresults.lowrankpart') ;
if ~exist(destDirlowrank,'dir')
    mkdir(destRoot,'PostMergePatchRPCAresults.lowrankpart') ;
end

destDirwithmask = fullfile(destRoot,'PostMergePatchRPCAresults.withmaskpart') ;
if ~exist(destDirwithmask,'dir')
    mkdir(destRoot,'PostMergePatchRPCAresults.withmaskpart') ;
end

[fileNames3, numImages3] = gzget_training_images( destRoot, 'afterfiltering') ;
% read every images
testImage = imread(fileNames3{1});

if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);

for fileIndex = 1:numImages3
    imPrepare= zeros(h,w,3, 'double');
    imPrepare(:,:,:) = 255;
    outputFileName  = sprintf('gztestlowrank%06d.jpg',fileIndex);
    outputFileNames= fullfile(destDirlowrank, outputFileName);     
    imwrite(uint8(imPrepare), outputFileNames);
    
    outputFileName  = sprintf('gztestwithmask%06d.jpg',fileIndex);
    outputFileNames= fullfile(destDirwithmask, outputFileName);     
    imwrite(uint8(imPrepare), outputFileNames);
end

%%change E_mat into binary matrix
for partnum = 1:6
    %clc ; clear all; close all ;
    MatrixERname  = sprintf('MatrixER%02d.mat',partnum*2-1);
    ERname = fullfile(destRoot2, MatrixERname); 
    %save(outputMatrixERname,'E_hatR') ;
    load(ERname);
    E_hatR1=E_hatR;
    clear E_hatR;
    
    MatrixERname  = sprintf('MatrixER%02d.mat',partnum*2);
    ERname = fullfile(destRoot2, MatrixERname); 
    %save(outputMatrixERname,'E_hatR') ;
    load(ERname);
    E_hatR2=E_hatR;
    clear E_hatR;
    E_hatR2part = [E_hatR1;E_hatR2];
    meanER = mean2(E_hatR2part);
    stdER  = std2(E_hatR2part);
    nRatio = 2;
    HighThreshR = meanER + nRatio*stdER;
    LowThreshR  = meanER - nRatio*stdER;
    CheckR1 = E_hatR2part < HighThreshR;
    CheckR2 = E_hatR2part > LowThreshR;
    CheckR = CheckR1 & CheckR2;
    CheckR = 1 - CheckR;
    
    MatrixEGname  = sprintf('MatrixEG%02d.mat',partnum*2-1);
    EGname = fullfile(destRoot2, MatrixEGname); 
    %save(outputMatrixEGname,'E_hatG') ;
    load(EGname);
    E_hatG1=E_hatG;
    
    MatrixEGname  = sprintf('MatrixEG%02d.mat',partnum*2);
    EGname = fullfile(destRoot2, MatrixEGname); 
    %save(outputMatrixEGname,'E_hatG') ;
    load(EGname);
    E_hatG2=E_hatG;
    clear E_hatG;
    E_hatG2part = [E_hatG1;E_hatG2];
    meanEG = mean2(E_hatG2part);
    stdEG  = std2(E_hatG2part);
    HighThreshG = meanEG + nRatio*stdEG;
    LowThreshG  = meanEG - nRatio*stdEG;
    CheckG1 = E_hatG2part < HighThreshG;
    CheckG2 = E_hatG2part > LowThreshG;
    CheckG = CheckG1 & CheckG2;
    CheckG = 1 - CheckG;

    
    MatrixEBname  = sprintf('MatrixEB%02d.mat',partnum*2-1);
    EBname = fullfile(destRoot2, MatrixEBname); 
    %save(outputMatrixEBname,'E_hatB') ;
    load(EBname);
    E_hatB1=E_hatB;
    
    MatrixEBname  = sprintf('MatrixEB%02d.mat',partnum*2);
    EBname = fullfile(destRoot2, MatrixEBname); 
    %save(outputMatrixEBname,'E_hatB') ;
    load(EBname);
    E_hatB2=E_hatB;
    clear E_hatB;
    
    E_hatB2part = [E_hatB1;E_hatB2];
    meanEB = mean2(E_hatB2part);
    stdEB  = std2(E_hatB2part);
    HighThreshB = meanEB + nRatio*stdEB;
    LowThreshB  = meanEB - nRatio*stdEB;
    CheckB1 = E_hatB2part < HighThreshB;
    CheckB2 = E_hatB2part > LowThreshB;
    CheckB = CheckB1 & CheckB2;
    CheckB = 1 - CheckB;
    CheckRGB=CheckR|CheckG|CheckB;

    [hRGB, wRGB] = size(CheckRGB);
    %     CheckRGB1=CheckRGB(1:hRGB/2,:);
    %     CheckRGB2=CheckRGB(hRGB/2+1:hRGB,:);
    EBinary=CheckRGB(1:hRGB/2,:);
    MatrixERGBname  = sprintf('MatrixEbinary%02d.mat',partnum*2-1);
    ERGBbinaryname = fullfile(destRoot2, MatrixERGBname); 
    save(ERGBbinaryname,'EBinary') ;
    
    EBinary=CheckRGB(hRGB/2+1:hRGB,:);
    MatrixERGBname  = sprintf('MatrixEbinary%02d.mat',partnum*2);
    ERGBbinaryname = fullfile(destRoot2, MatrixERGBname); 
    save(ERGBbinaryname,'EBinary');
    %destRoot2 = fullfile(destRoot,'afterfiltering') ;
end


for partnum = 1:12
    %clc ; clear all; close all ;
    [fileNameslowrank, numImageslowrank] = gzget_training_images( destRoot, 'PostMergePatchRPCAresults.lowrankpart') ;
    [fileNameswithmask, numImageswithmask] = gzget_training_images( destRoot, 'PostMergePatchRPCAresults.withmaskpart') ;
    [fileNamesOri, numImagesOri] = gzget_training_images( destRoot, 'afterfiltering') ;


    %destRoot2 = fullfile(destRoot,'afterfiltering') ;
    
%     MatrixERname  = sprintf('MatrixER%02d.mat',partnum);
%     ERname = fullfile(destRoot2, MatrixERname); 
%     %save(outputMatrixERname,'E_hatR') ;
%     load(ERname);
    
    MatrixARname  = sprintf('MatrixAR%02d.mat',partnum);
    ARname = fullfile(destRoot2, MatrixARname); 
    %save(outputMatrixARname,'A_hatR') ;
    load(ARname);
    
%     MatrixEGname  = sprintf('MatrixEG%02d.mat',partnum);
%     EGname = fullfile(destRoot2, MatrixEGname); 
%     %save(outputMatrixEGname,'E_hatG') ;
%     load(EGname);
    
    MatrixAGname  = sprintf('MatrixAG%02d.mat',partnum);
    AGname = fullfile(destRoot2, MatrixAGname); 
    %save(outputMatrixAGname,'A_hatG') ;
    load(AGname);
    

%     MatrixEBname  = sprintf('MatrixEB%02d.mat',partnum);
%     EBname = fullfile(destRoot2, MatrixEBname); 
%     %save(outputMatrixEBname,'E_hatB') ;
%     load(EBname);
    
    MatrixABname  = sprintf('MatrixAB%02d.mat',partnum);
    ABname = fullfile(destRoot2, MatrixABname); 
    %save(outputMatrixABname,'A_hatB') ;
    load(ABname);
    
    MatrixERGBname  = sprintf('MatrixEbinary%02d.mat',partnum);
    ERGBbinaryname = fullfile(destRoot2, MatrixERGBname); 
    %save(ERGBbinaryname,'EBinary') ;
    load(ERGBbinaryname);
    
%     FileNameMask  = sprintf('WithMask%02d.mat',partnum);
%     %load('WithMask.mat'); %save('WithMask.mat','WithMask') ;
%     load(FileNameMask);
% 
%     FileNameBK  = sprintf('Cleanbk%02d.mat',partnum);
%     % save('Cleanbk.mat','ClearBK') ;
%     load(FileNameBK);
    
    htemp=80;
    wtemp=80;
    
    for fileIndex = 1:numImageslowrank
        currentImageLowrank  = imread(fileNameslowrank{fileIndex});
        currentImageMask = imread(fileNameswithmask{fileIndex});
        
        OriImage = imread(fileNamesOri{fileIndex});
    
        BinaryMask = EBinary(:,fileIndex);
        BinaryMaskframe = reshape(BinaryMask,htemp,wtemp);
        
        tempBKR   = A_hatR(:,fileIndex);
        tempBKR = reshape(tempBKR,htemp,wtemp)*255;
            
        tempBKG   = A_hatG(:,fileIndex);
        tempBKG = reshape(tempBKG,htemp,wtemp)*255;
    
        tempBKB   = A_hatB(:,fileIndex);
        tempBKB = reshape(tempBKB,htemp,wtemp)*255;
        
%         gztest(:,:,1)=tempBKR;
%         gztest(:,:,2)=tempBKG;
%         gztest(:,:,3)=tempBKB;
%         imshow(gztest);
        
        if partnum==1
            currentImageLowrank(1:80,1:80,1) = tempBKR;
            currentImageLowrank(1:80,1:80,2) = tempBKG;
            currentImageLowrank(1:80,1:80,3) = tempBKB;
            
            Rori=double(OriImage(1:80,1:80,1));
            Gori=double(OriImage(1:80,1:80,2));
            Bori=double(OriImage(1:80,1:80,3));
            Rori=Rori.*(1.0 - BinaryMaskframe.*1.0) + BinaryMaskframe.*255.0;
            %Gori=Gori.*(1.0 - BinaryMaskframe.*1.0) + BinaryMaskframe.*255.0;
            Gori=Gori.*(1.0 - BinaryMaskframe.*1.0);
            Bori=Bori.*(1.0 - BinaryMaskframe.*1.0) + BinaryMaskframe.*255.0;
            
            currentImageMask(1:80,1:80,1) = Rori;
            currentImageMask(1:80,1:80,2) = Gori;
            currentImageMask(1:80,1:80,3) = Bori;
        end 
        if partnum==2
            currentImageLowrank(1:80,81:160,1) = tempBKR;
            currentImageLowrank(1:80,81:160,2) = tempBKG;
            currentImageLowrank(1:80,81:160,3) = tempBKB;
            
            Rori=double(OriImage(1:80,81:160,1));
            Gori=double(OriImage(1:80,81:160,2));
            Bori=double(OriImage(1:80,81:160,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(1:80,81:160,1) = Rori;
            currentImageMask(1:80,81:160,2) = Gori;
            currentImageMask(1:80,81:160,3) = Bori;
        end 
        if partnum==3
            currentImageLowrank(1:80,161:240,1) = tempBKR;
            currentImageLowrank(1:80,161:240,2) = tempBKG;
            currentImageLowrank(1:80,161:240,3) = tempBKB;
            
            Rori=double(OriImage(1:80,161:240,1));
            Gori=double(OriImage(1:80,161:240,2));
            Bori=double(OriImage(1:80,161:240,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(1:80,161:240,1) = Rori;
            currentImageMask(1:80,161:240,2) = Gori;
            currentImageMask(1:80,161:240,3) = Bori;
        end 
        if partnum==4
            currentImageLowrank(1:80,241:320,1) = tempBKR;
            currentImageLowrank(1:80,241:320,2) = tempBKG;
            currentImageLowrank(1:80,241:320,3) = tempBKB;
            
            Rori=double(OriImage(1:80,241:320,1));
            Gori=double(OriImage(1:80,241:320,2));
            Bori=double(OriImage(1:80,241:320,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(1:80,241:320,1) = Rori;
            currentImageMask(1:80,241:320,2) = Gori;
            currentImageMask(1:80,241:320,3) = Bori;
        end 
        
        if partnum==5
            currentImageLowrank(81:160,1:80,1) = tempBKR;
            currentImageLowrank(81:160,1:80,2) = tempBKG;
            currentImageLowrank(81:160,1:80,3) = tempBKB;
            
            Rori=double(OriImage(81:160,1:80,1));
            Gori=double(OriImage(81:160,1:80,2));
            Bori=double(OriImage(81:160,1:80,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(81:160,1:80,1) = Rori;
            currentImageMask(81:160,1:80,2) = Gori;
            currentImageMask(81:160,1:80,3) = Bori;
        end 
        if partnum==6
            currentImageLowrank(81:160,81:160,1) = tempBKR;
            currentImageLowrank(81:160,81:160,2) = tempBKG;
            currentImageLowrank(81:160,81:160,3) = tempBKB;
            
            Rori=double(OriImage(81:160,81:160,1));
            Gori=double(OriImage(81:160,81:160,2));
            Bori=double(OriImage(81:160,81:160,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(81:160,81:160,1) = Rori;
            currentImageMask(81:160,81:160,2) = Gori;
            currentImageMask(81:160,81:160,3) = Bori;
        end 
        if partnum==7
            currentImageLowrank(81:160,161:240,1) = tempBKR;
            currentImageLowrank(81:160,161:240,2) = tempBKG;
            currentImageLowrank(81:160,161:240,3) = tempBKB;
            
            Rori=double(OriImage(81:160,161:240,1));
            Gori=double(OriImage(81:160,161:240,2));
            Bori=double(OriImage(81:160,161:240,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe) ;
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(81:160,161:240,1) = Rori;
            currentImageMask(81:160,161:240,2) = Gori;
            currentImageMask(81:160,161:240,3) = Bori;
        end 
        if partnum==8
            currentImageLowrank(81:160,241:320,1) = tempBKR;
            currentImageLowrank(81:160,241:320,2) = tempBKG;
            currentImageLowrank(81:160,241:320,3) = tempBKB;
            
            Rori=double(OriImage(81:160,241:320,1));
            Gori=double(OriImage(81:160,241:320,2));
            Bori=double(OriImage(81:160,241:320,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe) ;
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(81:160,241:320,1) = Rori;
            currentImageMask(81:160,241:320,2) = Gori;
            currentImageMask(81:160,241:320,3) = Bori;
        end
        
        if partnum==9
            currentImageLowrank(161:240,1:80,1) = tempBKR;
            currentImageLowrank(161:240,1:80,2) = tempBKG;
            currentImageLowrank(161:240,1:80,3) = tempBKB;
            
            Rori=double(OriImage(161:240,1:80,1));
            Gori=double(OriImage(161:240,1:80,2));
            Bori=double(OriImage(161:240,1:80,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(161:240,1:80,1) = Rori;
            currentImageMask(161:240,1:80,2) = Gori;
            currentImageMask(161:240,1:80,3) = Bori;
        end 
        if partnum==10
            currentImageLowrank(161:240,81:160,1) = tempBKR;
            currentImageLowrank(161:240,81:160,2) = tempBKG;
            currentImageLowrank(161:240,81:160,3) = tempBKB;
            
            Rori=double(OriImage(161:240,81:160,1));
            Gori=double(OriImage(161:240,81:160,2));
            Bori=double(OriImage(161:240,81:160,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe) ;
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(161:240,81:160,1) = Rori;
            currentImageMask(161:240,81:160,2) = Gori;
            currentImageMask(161:240,81:160,3) = Bori;
        end 
        if partnum==11
            currentImageLowrank(161:240,161:240,1) = tempBKR;
            currentImageLowrank(161:240,161:240,2) = tempBKG;
            currentImageLowrank(161:240,161:240,3) = tempBKB;
            
            Rori=double(OriImage(161:240,161:240,1));
            Gori=double(OriImage(161:240,161:240,2));
            Bori=double(OriImage(161:240,161:240,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe) ;
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(161:240,161:240,1) = Rori;
            currentImageMask(161:240,161:240,2) = Gori;
            currentImageMask(161:240,161:240,3) = Bori;
        end 
        if partnum==12
            currentImageLowrank(161:240,241:320,1) = tempBKR;
            currentImageLowrank(161:240,241:320,2) = tempBKG;
            currentImageLowrank(161:240,241:320,3) = tempBKB;
            
            Rori=double(OriImage(161:240,241:320,1));
            Gori=double(OriImage(161:240,241:320,2));
            Bori=double(OriImage(161:240,241:320,3));
            Rori=Rori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            %Gori=Gori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            Gori=Gori.*(1 - BinaryMaskframe);
            Bori=Bori.*(1 - BinaryMaskframe) + BinaryMaskframe*255;
            
            currentImageMask(161:240,241:320,1) = Rori;
            currentImageMask(161:240,241:320,2) = Gori;
            currentImageMask(161:240,241:320,3) = Bori;
        end
        
        
        outputFileNameSparse  = sprintf('gztestlowrank%06d.jpg',fileIndex);
        outputFileNamesSparse= fullfile(destDirlowrank, outputFileNameSparse);
        imwrite(uint8(currentImageLowrank), outputFileNamesSparse);
        %imwrite(currentImageLowrank, outputFileNamesSparse);
        
        outputFileNameSparse  = sprintf('gztestwithmask%06d.jpg',fileIndex);
        outputFileNamesSparse= fullfile(destDirwithmask, outputFileNameSparse);
        imwrite(uint8(currentImageMask), outputFileNamesSparse);

%         outputFileName  = sprintf('gztestlowrank%06d.jpg',fileIndex);
%         outputFileNames= fullfile(destDir, outputFileName);
%         imwrite(uint8(imPrepare), outputFileNames);
%         outputFileName  = sprintf('gztestwithmask%06d.jpg',fileIndex);
%         outputFileNames= fullfile(destDir1, outputFileName);
%         imwrite(uint8(imPrepare), outputFileNames);

    end
    clear A_hatR;
    clear E_hatR;
    clear A_hatG;
    clear E_hatG;
    clear A_hatB;
    clear E_hatB;
end
