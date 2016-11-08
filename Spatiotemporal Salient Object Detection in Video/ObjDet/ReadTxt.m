
fin=fopen('2538-5_70133.tif.txt','r'); 
str = fgetl(fin);
Datas = strread(str, '%s');

GTBox = [str2double(Datas{1,1}),str2double(Datas{2,1}),str2double(Datas{3,1}),str2double(Datas{4,1})];

img=imread('F:\Action Recognition Data\UCF Sports\images\Diving-Side\001\2538-5_70133.jpg');
figure;imshow(img)
img(GTBox(2),GTBox(1):GTBox(1)+GTBox(3))= 255;
img(GTBox(2)+GTBox(4),GTBox(1):GTBox(1)+GTBox(3))= 255;
img(GTBox(2):GTBox(2)+GTBox(4),GTBox(1))= 255;
img(GTBox(2):GTBox(2)+GTBox(4),GTBox(1)+GTBox(3))= 255;
figure;imshow(img)
