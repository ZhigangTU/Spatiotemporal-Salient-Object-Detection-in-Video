%gz test purple color

imPrepare= zeros(240,320,3, 'double');
imPrepare(:,:,1) = 255;
imPrepare(:,:,2) = 0;
imPrepare(:,:,3) = 255;


% outputFileName  = sprintf('gztestlowrank%06d.jpg',fileIndex);
% outputFileNames= fullfile(destDirlowrank, outputFileName);
% imwrite(uint8(imPrepare), outputFileNames);
% 
% outputFileName  = sprintf('zzpurplecolor.jpg',fileIndex);
% outputFileNames= fullfile(destDirwithmask, outputFileName);
imshow(imPrepare);
imwrite(uint8(imPrepare), 'zzpurplecolor.jpg');