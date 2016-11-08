%% read every images
testImage = imread('SceneSampling000001.jpg');

I=double(testImage);
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);

Y = 0.299*R + 0.587*G + 0.114*B;
U = -0.147*R - 0.289*G + 0.436*B;
V = 0.615*R - 0.515*G - 0.100*B;
J=cat(3,Y,U,V);

% if isrgb(testImage)
%     grayImage = rgb2gray(testImage);
% end
% [h, w] = size(testImage);
% testImage = testImage(20:h,:,:);
% ReSizeIMG = imresize(testImage, [240 320]);

%imwrite(J, 'yuv000001.jpg');
imwrite(uint8(Y), 'Y000001.jpg');
imwrite(uint8(U), 'U000001.jpg');
imwrite(uint8(V), 'V000001.jpg');




