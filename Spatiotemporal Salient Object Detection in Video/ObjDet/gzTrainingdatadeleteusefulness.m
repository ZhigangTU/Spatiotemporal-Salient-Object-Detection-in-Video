%when there is too much training data, try to delete some usefulness

% clear
clc ; clear all; close all ;

addpath dataoftrainingROI ;
addpath dataoftraininggood ;
methodchoose =9001; 


currentPath = cd;
% input path
imagePath = fullfile(currentPath,'dataoftrainingROI') ;


%userName  = 'cam2new';
%userName  = 'cam4new';
userName  = 'cam1new';



destRoot = fullfile(currentPath,'dataoftraininggood');
destDir = fullfile(destRoot,userName) ;
if ~exist(destDir,'dir')
    mkdir(destRoot,userName) ;
end

% Get training images
[fileNames, numImages] = gzget_training_images( imagePath, userName) ;
usefulFlag = ones(numImages,1);
for fileIndex = 1 : numImages-1
    if usefulFlag(fileIndex)==1
        currentImage = imread(fileNames{fileIndex});
        if isrgb(currentImage)
            currentImage = rgb2gray(currentImage);
        end;
        
        %
        for fileIndexNext = fileIndex + 1 : numImages
            if usefulFlag(fileIndexNext)==1
                nextImage = imread(fileNames{fileIndexNext});
                if isrgb(nextImage)
                    nextImage = rgb2gray(nextImage);
                end;
                
                %
                imgDiff = abs(currentImage - nextImage);
                diffmaxvalue = max(max(imgDiff));
                
                if diffmaxvalue<8
                    usefulFlag(fileIndexNext)=0;
                end
                
            end
        end
        
    end

end

gaozhi = 100;
for fileIndex = 1:numImages

    if usefulFlag(fileIndex)==1
        currentImage = imread(fileNames{fileIndex});

        outputname  = sprintf('usefultraining%d.jpg',fileIndex);
        outputnames = fullfile(destDir, outputname);
        
        imwrite(currentImage, outputnames);
        %imwrite(uint8(BkImage), outputFileNamesSparse);
    end
    
end