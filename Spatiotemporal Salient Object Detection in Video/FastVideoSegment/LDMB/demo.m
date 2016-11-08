% demo for detecting motion boundaries with different cues, please see readme.txt first
%
% Piotr's Toolbox and Classic+NL flow code
% addpath(genpath(<path_to_piotr_toolbox>));
% addpath(genpath(<path_to_classicnl_flow>));


image = imread('im1.jpg');
next_image = imread('im2.jpg');
prev_image = imread('im0.jpg');

% forward flow computation
load('flow.mat'); % flow = estimate_flow_interface(image, next_image, 'classic+nl-fast');

load('backward_flow.mat'); % backward_flow = estimate_flow_interface(image, prev_image, 'classic+nl-fast');

% using color only
boundaries_Color = detect_motionboundaries(image); % not accurate

% using color and optical flow
boundaries_ColorFlow = detect_motionboundaries(image, flow);

% using color, optical flow and warping error
boundaries_ColorFlowWarping = detect_motionboundaries(image, flow, next_image); % recommanded

% using color, optical flow, warping error and backward information
boundaries_ColorFlowWarpingBackward = detect_motionboundaries(image, flow, next_image, backward_flow, prev_image); % recommanded but slower since it requires backward flow

% display
hold on;
subplot(221);
imshow(boundaries_Color);
title('Color');
subplot(222);
imshow(boundaries_ColorFlow);
title('Color+Flow');
subplot(223);
imshow(boundaries_ColorFlowWarping);
title('Color+Flow+Warping');
subplot(224);
imshow(boundaries_ColorFlowWarpingBackward);
title('Color+Flow+Warping+Backward');
hold off;
