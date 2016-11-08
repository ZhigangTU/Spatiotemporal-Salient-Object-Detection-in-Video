% demo for detecting motion boundaries with different cues, please see readme.txt first
%
% Piotr's Toolbox and Classic+NL flow code
% addpath(genpath(<path_to_piotr_toolbox>));
% addpath(genpath(<path_to_classicnl_flow>));

outPath = 'E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images\Walk-Front\003\LDMB';
image = imread('E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images\Walk-Front\003\RF1-10578_70289.jpg');
next_image = imread('E:\Action Recognition\Pose-CNN\Data\UCF-Sports\images\Walk-Front\003\RF1-10578_70290.jpg');

% % forward flow computation
% load('flow.mat'); % flow = estimate_flow_interface(image, next_image, 'classic+nl-fast');
load('E:\Action Recognition\Pose-CNN\cache\OFUCFSport\Walk-Front\003\flowu5.mat');
load('E:\Action Recognition\Pose-CNN\cache\OFUCFSport\Walk-Front\003\flowv5.mat');
flow = cat(3,u,v);

% %backward flow
% prev_image = imread('im0.jpg');
% load('backward_flow.mat'); % backward_flow = estimate_flow_interface(image, prev_image, 'classic+nl-fast');

% using color only
boundaries_Color = detect_motionboundaries(image); % not accurate
figure; imshow(boundaries_Color)
figure; imshow(boundaries_Color>0.05)
figure; imshow(boundaries_Color>0.1)
outObjectPath = sprintf('%s/boundaries_Color-Fram005.jpg',outPath);    
imwrite(boundaries_Color, outObjectPath);
boundaries_ColorT005 = boundaries_Color>0.05;
outObjectPath = sprintf('%s/boundaries_ColorT005-Fram005.jpg',outPath);    
imwrite(boundaries_ColorT005, outObjectPath);
boundaries_ColorT01 = boundaries_Color>0.1;
outObjectPath = sprintf('%s/boundaries_ColorT01-Fram005.jpg',outPath);    
imwrite(boundaries_ColorT01, outObjectPath);

% using color and optical flow
boundaries_ColorFlow = detect_motionboundaries(image, flow);
figure; imshow(boundaries_ColorFlow)
figure; imshow(boundaries_ColorFlow>0.05)
figure; imshow(boundaries_ColorFlow>0.1)
outObjectPath = sprintf('%s/boundaries_ColorFlow-Fram005.jpg',outPath);
imwrite(boundaries_ColorFlow, outObjectPath);
boundaries_ColorFlowT005 = boundaries_ColorFlow>0.05;
outObjectPath = sprintf('%s/boundaries_ColorFlowT005-Fram005.jpg',outPath);    
imwrite(boundaries_ColorFlowT005, outObjectPath);
boundaries_ColorFlowT01 = boundaries_ColorFlow>0.1;
outObjectPath = sprintf('%s/boundaries_ColorFlowT01-Fram005.jpg',outPath);    
imwrite(boundaries_ColorFlowT01, outObjectPath);

% using color, optical flow and warping error
boundaries_ColorFlowWarping = detect_motionboundaries(image, flow, next_image); % recommanded
figure; imshow(boundaries_ColorFlowWarping)
figure; imshow(boundaries_ColorFlowWarping>0.05)
figure; imshow(boundaries_ColorFlowWarping>0.1)
outObjectPath = sprintf('%s/boundaries_ColorFlowWarping-Fram005.jpg',outPath);
imwrite(boundaries_ColorFlowWarping, outObjectPath);
boundaries_ColorFlowWarpingT005 = boundaries_ColorFlowWarping>0.05;
outObjectPath = sprintf('%s/boundaries_ColorFlowWarpingT005-Fram005.jpg',outPath);    
imwrite(boundaries_ColorFlowWarpingT005, outObjectPath);
boundaries_ColorFlowWarpingT01 = boundaries_ColorFlowWarping>0.1;
outObjectPath = sprintf('%s/boundaries_ColorFlowWarpingT01-Fram005.jpg',outPath);    
imwrite(boundaries_ColorFlowWarpingT01, outObjectPath);

% % using color, optical flow, warping error and backward information
% boundaries_ColorFlowWarpingBackward = detect_motionboundaries(image, flow, next_image, backward_flow, prev_image); % recommanded but slower since it requires backward flow

% % display
% hold on;
% subplot(221);
% imshow(boundaries_Color);
% title('Color');
% subplot(222);
% imshow(boundaries_ColorFlow);
% title('Color+Flow');
% subplot(223);
% imshow(boundaries_ColorFlowWarping);
% title('Color+Flow+Warping');
% subplot(224);
% imshow(boundaries_ColorFlowWarpingBackward);
% title('Color+Flow+Warping+Backward');
% hold off;
