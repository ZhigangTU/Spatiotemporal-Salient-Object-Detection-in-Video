% Gao Zhi rpca for color sequence, to obatin background
% robust batch image alignment example

% clear
clc ; clear all; close all ;

OriImage = imread('Frame0003201.jpg');
OriR = double(OriImage(:,:,1));
OriG = double(OriImage(:,:,2));
OriB = double(OriImage(:,:,3));
  
BKImage  = imread('background1.jpg');
BKR = double(BKImage(:,:,1));
BKG = double(BKImage(:,:,2));
BKB = double(BKImage(:,:,3));

DiffR=abs(OriR-BKR)/255;
DiffG=abs(OriG-BKG)/255;
DiffB=abs(OriB-BKB)/255;
Maxdiff=max(max(DiffR,DiffG),DiffB);

imshow(Maxdiff);









