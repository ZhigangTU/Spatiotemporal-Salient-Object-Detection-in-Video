
% h=100;
% w=255*5;
% 
% dataR = zeros([h, w], 'double');
% dataG = zeros([h, w], 'double');
% dataB = zeros([h, w], 'double');
% 
% for fileIndex = 1:1:255
%   dataR(:,(fileIndex-1)*5+1 : (fileIndex-1)*5+5)= fileIndex;
% end
%   
% BkImage(:,:,1)=dataR;
% BkImage(:,:,2)=dataG;
% BkImage(:,:,3)=dataB;
% 
% % utputFileNameSparse  = sprintf('background%d.jpg',fileIndex);
% % outputFileNamesSparse= fullfile(destDir, outputFileNameSparse);
% %imwrite(im2, outputFileNamesSparse);
% %imwrite(uint8(im2), outputFileNamesSparse);
% % imwrite(uint8(BkImage), outputFileNamesSparse);
% imwrite(uint8(BkImage), 'generatered.bmp');


% IMin0 = imread(strcat([pathForImages,imgLibrary{img},'.png']));
IMin0 = imread('gzoutputmask0080.PNG');
dataR0=IMin0(:,:,1);

[row1, col1] = find( dataR0~= 0 );
% MatNonzero= dataR0(row1, col1);
num = size(row1, 1);            % ????????
MinValue=255;
Valuevector=[];
for i = 1:num
    Valuevector=[Valuevector dataR0(row1(i), col1(i))];
    if dataR0(row1(i), col1(i))<MinValue
        MinValue=dataR0(row1(i), col1(i));
    end
end
mm=min(Valuevector);



maskR = dataR0>10;
maskRInverse = 1 - maskR;

%maskRInverse = double (maskRInverse);
% dataG0=IMin0(:,:,2);
% dataB0=IMin0(:,:,3);

IMin1 = imread('frame3_0080.bmp');

% testImage = imread(fileNames{1});

if isrgb(IMin1)
    IMin1gray = rgb2gray(IMin1);
end
% [h, w] = size(IMin1gray);
dataOri = double (IMin1gray);
dataOriSelect = dataOri.* maskR;

[row2, col2] = find( dataOriSelect~= 0 );
num2 = size(row2, 1);            % ????????
MinValue2=255;
Valuevector2=[];
for i = 1:num2
    Valuevector2=[Valuevector2 dataOriSelect(row2(i), col2(i))];
    if dataOriSelect(row2(i), col2(i))<MinValue2
        MinValue2=dataOriSelect(row2(i), col2(i));
    end
end
mm2=min(Valuevector2);
mm3=max(Valuevector2);

fRatio = (double (mm))/(double (mm3));
% dataOriSelectReset = (dataOri.* maskR)*fRatio;
dataOriSelectReset = dataOriSelect*fRatio;

dataR0 = double (dataR0);
IMinnew(:,:,1)=dataOri.*maskRInverse + dataR0;
IMinnew(:,:,2)=dataOri.*maskRInverse + dataOriSelectReset;
IMinnew(:,:,3)=dataOri.*maskRInverse + dataOriSelectReset;


% dataR1=IMin1(:,:,1);
% dataG1=double (IMin1(:,:,2));
% dataB1=double (IMin1(:,:,3));
% 
% dataR1=double (dataR1);
% dataR0=double (dataR0);
% dataR1New = dataR1.*maskRInverse+dataR0;
% % dataR1New = dataR1+dataR0;
% 
% IMinnew(:,:,1)=dataR1New;
% IMinnew(:,:,2)=dataG1;
% IMinnew(:,:,3)=dataB1;
imwrite(IMinnew/255, 'salient.bmp');
gaozhi=0;
