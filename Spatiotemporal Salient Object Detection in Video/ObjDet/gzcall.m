% Figure 5 in the paper
% robust batch image alignment example

% clear all;
% close all ;

% addpath
addpath computeopticflow ;
addpath data ;
addpath results ;

%% define images' path
currentPath = cd;

% input path
imagePath = fullfile(currentPath,'data') ;
pointPath = fullfile(currentPath,'data') ; % path to files containing initial feature coordinates
userName = 'carheavyrain';

% output path
destRoot = fullfile(currentPath,'results') ;
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

%% Get training images
[fileNames, numImages] = gzget_training_images(imagePath, userName) ;

%% read every images
testImage = imread(fileNames{1});
% double(imread(fileNames{fileIndex}));
[h, w] = size(testImage);
   
%for color image
w=w/3;
%norm_matrix = zeros(h*w,numImages-1);

for fileIndex = 1 : numImages-1
    
    currentImage = imread(fileNames{fileIndex});%double(imread(fileNames{fileIndex}));
    nextImage    = imread(fileNames{fileIndex+1});%double(imread(fileNames{fileIndex+1}));
    
    uv = estimate_flow_interface(currentImage, nextImage,'classic+nl-fast');
    u = uv(:,:,1);
    v = uv(:,:,2);
    
    norm=zeros(h,w);
    
   for j = 1:w 
       for i = 1:h 
           norm(i,j) = sqrt(u(i,j) ^ 2 + v(i,j) ^ 2);
       end
   end
   
   img = uint8(flowToColor(uv));
   
   NormName = sprintf('norm%d.mat',fileIndex);
   NormNames = fullfile(destDir, NormName); 
   save(NormNames,'norm') ;
   
   UName = sprintf('flowu%d.mat',fileIndex);
   UNames = fullfile(destDir, UName); 
   save(UNames,'u') ;
   
   VName = sprintf('flowv%d.mat',fileIndex);
   VNames = fullfile(destDir, VName); 
   save(VNames,'v') ;
   
   outputFileName = sprintf('tree%d.jpg',fileIndex);
   outputFileNames = fullfile(destDir, outputFileName); 
   imwrite(img, outputFileNames );
      
   %%%%%%norm_matrix(:,fileIndex) = norm;   
end

