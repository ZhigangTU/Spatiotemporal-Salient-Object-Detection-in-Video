% Gao Zhi rpca for color sequence, to obatin background
% robust batch image alignment example

% clear
clc ; clear all; close all ;

% addpath
addpath datacompare1 ;
addpath datacompare2 ;

addpath results ;
%%
methodchoose =1, %original RPCA code
%% define images' path
currentPath = cd;

% input path
imagePath1 = fullfile(currentPath,'datacompare1') ;
imagePath2 = fullfile(currentPath,'datacompare2') ;

%userName = 'manyBirds';
%userName = '5ready';
%userName = '2ready';
userName = 'Frontier02on20130118outputbmpchange';

% output path
if methodchoose==1
    destRoot = fullfile(currentPath,'datacompareresults') ;
end   

destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNames1, numImages1] = gzget_training_images( imagePath1, userName) ;
[fileNames2, numImages2] = gzget_training_images( imagePath2, userName) ;

if (numImages1~=numImages2)
    disp('the image numbers of these two folders are not equal! return!') ;
    return;
end 


outputTXT= fullfile(destDir, 'compareresults.txt');
%outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);

fid = fopen(outputTXT, 'wt');



for fileIndex = 1:numImages1
  Image1 = imread(fileNames1{fileIndex});
  Image1R = double (Image1(:,:,1));
  [h1, w1] = size(Image1R);
  
  Image1B = Image1(:,:,3);
  
  Image2 = imread(fileNames2{fileIndex});
  Image2R = double (Image2(:,:,1));
  [h2, w2] = size(Image2R);
  
  if (h1~=h2 | w1~=w2)
      disp('the dims of these two images do not match! return!') ;
      return;
  end

%   imshow(Image1R);
%   imshow(Image1B);
%   
%  
%   imshow(Image2R);
%   imshow(Image2R);

  diffimg=Image1R-Image2R;
%   imshow(Image1R);
%   imshow(Image2R);
%   imshow(diffimg);
  
  Absdiff=abs(diffimg);
  BinAbs=Absdiff>2;
  
  Number=sum(BinAbs(:)); %????
  fprintf(fid, '%06d %d\n', fileIndex, Number);
  
  Abssave=BinAbs*255;
  %  imshow(diffimg);
  %  imshow(diffimg);


  outputFileNameSparse  = sprintf('compare%06d.jpg',fileIndex);
  outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);
  %imwrite(im2, outputFileNamesSparse);
  %imwrite(uint8(im2), outputFileNamesSparse);
  imwrite(uint8(Abssave), outputFileNamesSparse);
end

fclose(fid);

% %% read every images
% testImage = imread(fileNames{1});
% 
% if isrgb(testImage)
%     testImage = rgb2gray(testImage);
% end
% [h, w] = size(testImage);
% 
% dataR = zeros([h*w, numImages], 'double');
% dataG = zeros([h*w, numImages], 'double');
% dataB = zeros([h*w, numImages], 'double');
% for fileIndex = 1:numImages
%   currentImage = imread(fileNames{fileIndex});
%   
%   currentImageR = currentImage(:,:,1);
%   currentImageG = currentImage(:,:,2);
%   currentImageB = currentImage(:,:,3);
%   
% %   if isrgb(currentImage)
% %     currentImage = rgb2gray(currentImage);
% %   end
%   
%   imR = reshape(currentImageR,w*h,1);
%   dataR(:, fileIndex) = imR(:,1);
%   
%   imG = reshape(currentImageG,w*h,1);
%   dataG(:, fileIndex) = imG(:,1);
%   
%   imB = reshape(currentImageB,w*h,1);
%   dataB(:, fileIndex) = imB(:,1);
% end
% 
% dataR=dataR/255.0;
% dataG=dataG/255.0;
% dataB=dataB/255.0;
% 
% if methodchoose==1
%     [A_hatR E_hatR iterR] = inexact_alm_rpca(dataR);%original ALM method
%     clear E_hatR;
%     clear iterR;
%     clear dataR;
%     
%     [A_hatG E_hatG iterG] = inexact_alm_rpca(dataG);%original ALM method
%     clear E_hatG;
%     clear iterG;
%     clear dataG;
%     
%     [A_hatB E_hatB iterB] = inexact_alm_rpca(dataB);%original ALM method
%     clear E_hatB;
%     clear iterB;
%     clear dataB;
%     
%     %[A_hat E_hat iter] = inexact_alm_rpca(data, 1 / sqrt(160*120));%original ALM method
% end   
% 
% %%save the data as images (low rank and sparse matrix)
% for fileIndex = 1:numImages
% %     immask = Mask_matrix(:,fileIndex);
% %     immask = reshape(immask,h,w);
% %     immask = immask * 255;
% %     %not display now
% %     %figure;
% %     %imshow(uint8(im1));
% %     outputMask  = sprintf('mask%d.jpg',fileIndex);
% %     outputMasks = fullfile(destDir, outputMask); 
% %     imwrite(immask, outputMasks);  
%     imR = A_hatR(:,fileIndex);
%     imR = reshape(imR,h,w);
%     imR = imR * 255;
%     
%     imG = A_hatG(:,fileIndex);
%     imG = reshape(imG,h,w);
%     imG = imG * 255;
%     
%     imB = A_hatB(:,fileIndex);
%     imB = reshape(imB,h,w);
%     imB = imB * 255;
%     
%     BkImage(:,:,1)=imR;
%     BkImage(:,:,2)=imG;
%     BkImage(:,:,3)=imB;
%     
%     outputFileNameSparse  = sprintf('background%d.jpg',fileIndex);
%     outputFileNamesSparse= fullfile(destDir, outputFileNameSparse); 
%     %imwrite(im2, outputFileNamesSparse);
%     %imwrite(uint8(im2), outputFileNamesSparse);
%     imwrite(uint8(BkImage), outputFileNamesSparse);
% end
