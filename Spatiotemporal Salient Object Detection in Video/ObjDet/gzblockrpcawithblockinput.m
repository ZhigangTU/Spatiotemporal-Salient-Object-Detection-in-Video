% robust batch image alignment example

% clear
clc ;
clear all;
close all ;

% addpath
addpath data ;
addpath results ;
addpath oriimagedata ;
addpath blockdata ;

%%
%methodchoose =1, %original RPCA code
%methodchoose =5, 
%             =2, try to enforce Block constraint
%methodchoose  =3;% basic mix 2 and 1 norm regularizer
%methodchoose  =4,
%methodchoose   =6,% block 2 1 norm 
%methodchoose   =8,
%fix block with known size and position
%methodchoose   = 11,
%methodchoose   = 12,
%methodchoose   = 13,

methodchoose   = 100,
%% define images' path
currentPath = cd;

% input path
%imagePath = fullfile(currentPath,'data') ;
%pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates

imagePath = fullfile(currentPath,'oriimagedata') ;
blockPath = fullfile(currentPath,'blockdata') ; % path to files containing initial feature coordinates

%userName = 'simulateblock6by10up120' ;
userName = 'smtreeman97320' ;

% output path
if methodchoose==11
    destRoot = fullfile(currentPath,'gzfixblockknowpositionandsizejustprpcesssuchblocks') ;
end 

if methodchoose==12
    destRoot = fullfile(currentPath,'gzfixblockknowpositionandsizeprocessallthensuchblocks') ;
end 

if methodchoose==13
    destRoot = fullfile(currentPath,'gzjustregularprocess') ;
end 

if methodchoose==100
    destRoot = fullfile(currentPath,'block rpca with already estimated block') ;
end 

if methodchoose==2
    destRoot = fullfile(currentPath,'gzblockrpcaresults') ;
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

data = zeros([h*w, numImages], 'double');
for fileIndex = 1:numImages
  currentImage = imread(fileNames{fileIndex});
  
  if isrgb(currentImage)
    currentImage = rgb2gray(currentImage);
  end
  
  im = reshape(currentImage,w*h,1);
  data(:, fileIndex) = im(:,1);
end

data=data/255.0;

if methodchoose==8
    [A_hat E_hat iter] = gzfixblock2and1norm_inexact_alm_rpca(data);
    %[A_hat E_hat iter] = gzblock2and1norm_inexact_alm_rpca(data);%
    %[Mask_matrix E_blockmask]=gzAnalyzeEmatrix(E_hat);
    %[A_hat E_hat iter] = gzAdaptiveThreshblock2and1norm_inexact_alm_rpca(data, E_blockmask);%
end 

if methodchoose==11
    [A_hat E_hat iter] = gzKnownblockparts2and1norm_inexact_alm_rpca(data);
end 

if methodchoose==12
    [A_hat E_hat iter] = gzKnownblockall2and1norm_inexact_alm_rpca(data);
end 

if methodchoose==100
    [A_hat E_hat iter] = gzBlockrpcaWithestimatedBlockinfo(data);
end 

if methodchoose==13
    [A_hat E_hat iter] = gzJustregularlockall2and1norm_inexact_alm_rpca(data);
end 

AANames = fullfile(destDir, 'Amatrix.mat'); 
save(AANames,'A_hat') ;
EENames = fullfile(destDir, 'Ematrix.mat'); 
save(EENames,'E_hat') ;

%%save the data as images (low rank and sparse matrix)
for fileIndex = 1:numImages
    %immask = Mask_matrix(:,fileIndex);
    %immask = reshape(immask,h,w);
    %immask = immask * 255;
    %outputMask  = sprintf('mask%d.jpg',fileIndex);
    %outputMasks = fullfile(destDir, outputMask); 
    %imwrite(immask, outputMasks);
    
    im1 = E_hat(:,fileIndex);
    for j = 1:h*w,
        %if (abs(im1(j,1))) < 0.02,
	     %    im1(j,1) = 0;
         %else im1(j,1)=255;
        %end
        
        if im1(j,1) < 0,
	      im1(j,1) = -1*im1(j,1);
        end
    end
    im1 = reshape(im1,h,w);
    %not display now
    %figure;
    %imshow(uint8(im1));
    
    outputFileNameLowrank  = sprintf('Lowrank%05d.jpg',fileIndex);
    outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
    imwrite(im1, outputFileNamesLowrank);
    
    im2 = A_hat(:,fileIndex);
    im2 = reshape(im2,h,w);
    im2 = im2 * 255;
    %not display now
    %figure;
    %imshow(uint8(im2));
    
    outputFileNameSparse  = sprintf('Sparse%05d.jpg',fileIndex);
    outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
    %imwrite(im2, outputFileNamesSparse);
    imwrite(uint8(im2), outputFileNamesSparse);
end

% 2nd ronund decomposition  %data=A_hat;