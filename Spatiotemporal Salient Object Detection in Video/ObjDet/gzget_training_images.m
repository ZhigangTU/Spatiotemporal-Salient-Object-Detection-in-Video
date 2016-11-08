function [fileNames, numImages] = gzget_training_images(rootPath, trainingDatabaseName)
% (imagePath, pointPath, userName)

fileNames = {};
imageIndex = 0;
userDirectoryContents = list_image_files(fullfile(rootPath, trainingDatabaseName));

 for fileIndex = 1:length(userDirectoryContents),
        imageName = userDirectoryContents{fileIndex};
%         disp(['Using image file ' imageName '...']);

        imageIndex = imageIndex+1;

        imageFileName = fullfile(rootPath, trainingDatabaseName, imageName);
        fileNames{imageIndex} = imageFileName;
 end
 
numImages = length(userDirectoryContents);