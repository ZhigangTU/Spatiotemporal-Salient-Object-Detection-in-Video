function [ boundaries ] = detect_motionboundaries( image, flow, next_image, backward_flow, prev_image, ms, nms, sharpen )
% Detect motion boundaries in an image
%
% For an introductory tutorial, please see demo.m
%
% INPUTS
% image                     image from which motion boundaries are extracted (H x W x 3 array)
% flow             []       forward flow to the next frame  (H x W x 2 array or []. if [], extract boundaries based on RGB only)
% next_image       []       next image in the video (H x W x 3 array or []. if [], will extract boundaries based on RGB and flow only)
% backward_flow    []       backward flow to the previous frame (H x W x 2 array or []. if [], will extract boundaries based on RGB, flow and image warping error only)
% prev_image       []       previous image in the video (H x W x 3 array if backward_flow is not empty, otherwise, empty array)
% ms               1        if true perform detection at multiple scale
% nms              0        if true apply non-maximum suppression to edges
% sharpen          2        sharpening amount (must be smaller than 2)
%
% OUPUTS
% boundaries            detected motion boundaries
%
% Code originally written by Philippe Weinzaepfel, 2015, 
% for motion boundaries detection
% Licensed under the MSR-LA Full Rights License [see license.txt]

if nargin<2, flow = []; end;
if nargin<3, next_image = []; end;
if nargin<4, backward_flow = []; end;
if nargin<5, prev_image = []; end;
if nargin<6, ms = 1; end;
if nargin<7, nms = 0; end;
if nargin<8, sharpen = 2; end;

% prepare data
inputs = 'Color';
I = image;
if ~isempty(flow),
    inputs = 'Color+Flow';
    I = cat(3, single(I)/255.0, single(flow)/10.0);
    if ~isempty(next_image),
        inputs = 'Color+Flow+Warping';
        [lab_error, hog_error] = image_warping_error(single(image), single(next_image), single(flow));
        I = cat(3, I, lab_error, hog_error);
        if ~isempty(backward_flow),
            inputs = 'Color+Flow+Warping+Backward';
            [lab_error, hog_error] = image_warping_error(single(image), single(prev_image), single(backward_flow));
            I = cat(3, I, single(backward_flow)/10.0, lab_error, hog_error);
        end;
    end;
end;

% load model
load( sprintf('MotionBoundariesModels/forest/model_SintelClean_Classic+NL-Fast_%s.mat', inputs));

model.opts.multiscale = ms;
model.opts.nms = nms;
model.opts.sharpen = sharpen;

% run motion boundaries detection
boundaries = edgesDetect(I, model);


end

