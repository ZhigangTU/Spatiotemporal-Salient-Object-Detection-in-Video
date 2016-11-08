% Figure 5 in the paper
% robust batch image alignment example

% clear
clc ;
clear all;
close all ;

% addpath
addpath data ;
addpath results ;
%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;
pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates
userName = 'smtreeman100320' ;

% output path
destRoot = fullfile(currentPath,'iterativeALM') ;


%% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;
%% read every images
testImage = imread(fileNames{1});
if isrgb(testImage)
    testImage = rgb2gray(testImage);
end
[h, w] = size(testImage);
   
%for color image
%w=w/3;

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

for nIteration = 1:20
    outputIteration = sprintf('of No %d iteration',nIteration);
    strName=fullfile(userName,outputIteration) ; 
    destDir = fullfile(destRoot,strName) ;
    if ~exist(destDir,'dir')
        mkdir(destRoot,strName) ;
    end
    mylamda=1*0.9.^(nIteration-1)/ sqrt(h*w);
    %[A_hat E_hat iter] = gzModifyblock_inexact_alm_rpca(data);%try to enforce block constraint,improve the previous case
    [A_hat E_hat iter] = gzblock_inexact_alm_rpca(data,mylamda);%try to enforce block constraint,improve the previous case
    
    AANames = fullfile(destDir, 'Amatrix.mat'); 
    save(AANames,'A_hat') ;
    EENames = fullfile(destDir, 'Ematrix.mat'); 
    save(EENames,'E_hat') ;
    
    %%save the data as images (low rank and sparse matrix)
    for fileIndex = 1:numImages
    
      im1 = E_hat(:,fileIndex);
      for j = 1:h*w,
          if (abs(im1(j,1))) < 0.02,
	          im1(j,1) = 0;
          else im1(j,1)=255;
          end
        
          %if im1(j,1) ~= 0,
	      %   im1(j,1) = 255;
          %end
      end
      im1 = reshape(im1,h,w);
      %not display now
      %figure;
      %imshow(uint8(im1));
    
      outputFileNameLowrank  = sprintf('Lowrank%d.jpg',fileIndex);
      outputFileNamesLowrank = fullfile(destDir, outputFileNameLowrank); 
      imwrite(im1, outputFileNamesLowrank);
    
    
      im2 = A_hat(:,fileIndex);
      im2 = reshape(im2,h,w);
      im2 = im2 * 255;
      %not display now
      %figure;
      %imshow(uint8(im2));
     
      outputFileNameSparse  = sprintf('Sparse%d.jpg',fileIndex);
      outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
      %imwrite(im2, outputFileNamesSparse);
      imwrite(uint8(im2), outputFileNamesSparse);
   end

   data=A_hat;
end   

%[A_hat E_hat iter] = inexact_alm_rpca(data);%original ALM method
%[A_hat E_hat iter] = gzblock_inexact_alm_rpca(data);%try to
%[A_hat E_hat iter] = gzModifyblock_inexact_alm_rpca(data);%try to 
%[A_hat E_hat iter] = gzbasic2and1_inexact_alm_rpca(data);%try to 
%[A_hat E_hat iter] = gzsimpleblock_inexact_alm_rpca(data);%try to





