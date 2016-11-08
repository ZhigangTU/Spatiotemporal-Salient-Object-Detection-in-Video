% Function to create segmentation visualisations
%
%    Copyright (C) 2013  Anestis Papazoglou
%
%    You can redistribute and/or modify this software for non-commercial use
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    For commercial use, contact the author for licensing options.
%
%    Contact: a.papazoglou@sms.ed.ac.uk

function createSegmentationVideo(options, params, data, mode)
    folder = fullfile(options.outfolder, 'visuals', 'segmentations');
    if(~exist(folder, 'dir'))
        mkdir(folder);
    end
    
    if(isfield(data, 'id') && ~isempty(data.id))
        range = data.id;
    else
        range = 1;
    end
    
    % Check that mode is set, otherwise set to default
    if(~exist('mode', 'var') || isempty(mode))
        mode = 'BoundsOverlay';
    end
    
    % Create video object
    if(isfield(params, 'name') && ~isempty(params.name))
        filename = fullfile(folder, sprintf('%s.avi', params.name));
    else
        filename = fullfile(folder, ...
            sprintf('segmentation-%s-%03d.avi', mode, range));
    end
    writerObject = VideoWriter(filename);
    if(isfield(params, 'frameRate') && ~isempty(params.frameRate))
        writerObject.FrameRate = params.frameRate;
    else
        writerObject.FrameRate = 15;
    end
    if(isfield(params, 'quality') && ~isempty(params.quality))
        writerObject.Quality = params.quality;
    else
        writerObject.Quality = 50;
    end
    open(writerObject);
    
    framesdata = data.imgs;
    segmentation = data.segmentation;
    
    if(strcmp(mode, 'BoundsOverlay'))
        for(index = 1: length(segmentation))
            frameid = options.ranges(range) + index - 1;
            frame = readFrame(options, frameid);

            segImg = overlaySegmentationBoundary(frame, ...
                segmentation{ index });
            fullframe = segImg;
            writeVideo(writerObject, fullframe);
        end
    elseif(strcmp(mode, 'BoundsInOutOverlay'))
        inMaps = data.inMaps;
        for(index = 1: length(segmentation) - 1)
            frameid = options.ranges(range) + index - 1;
            frame = readFrame(options, frameid);
            
            segImg = overlaySegmentationBoundary(frame, ...
                segmentation{ index });
            inoutImg = overlayInMap(frame, inMaps{ index });
            fullframe = [ frame, inoutImg, segImg ];
            writeVideo(writerObject, fullframe);
        end
    elseif(strcmp(mode, 'BoundsFlowInOutOverlay'))
        inMaps = data.inMaps;
        flow = data.flow;
        cutMask = false(size(inMaps{ 1 }));
        cutMask(1: 20, :) = true;
        cutMask(end - 20: end, :) = true;
        cutMask(:, 1: 20) = true;
        cutMask(:, end - 20: end) = true;
        for(index = 1: length(segmentation) - 1)
            frameid = options.ranges(range) + index - 1;
            frame = readFrame(options, frameid);
            
            mbounds = getProbabilityEdge(flow{ index }, 3) > 0.5;
            mbounds = cleanBoundaries(mbounds, cutMask);
            
            segImg = overlaySegmentationBoundary(frame, ...
                segmentation{ index });
            inoutImg = overlayInMap(frame, inMaps{ index });
            inoutImg = overlayMotionBoundary(inoutImg, mbounds);
            flowImg = flowToColor(flow{ index }, 15);
            fullframe = [ frame, flowImg; inoutImg, segImg ];
            fullframe = imresize(fullframe, 720 / size(fullframe, 1));
            writeVideo(writerObject, fullframe);
        end
        
    elseif(strcmp(mode, 'MaskOverlay'))
        for(index = 1: length(segmentation))
            frameid = options.ranges(range) + index - 1;
            frame = readFrame(options, frameid);
            
            segImg = overlaySegmentationMask(frame, ...
                segmentation{ index });
            fullframe = [ frame, segImg ];
            writeVideo(writerObject, fullframe);
        end
        
    elseif(strcmp(mode, 'ShowProcess'))
        locationProbability = data.locationProbability;
        appearanceProbability = data.appearanceProbability;
        unaryPotential = data.unaryPotential;
        for index = 1: length(segmentation) - 1  %TODO remove -1
            frame = framesdata{index};   % readFrame(options, index);
            segImg = overlaySegmentationMask(frame, segmentation{ index });
            
%             locationImg = gray2rgb(255 * locationProbability{ index });
%             appearanceImg = gray2rgb(255 * appearanceProbability{ index });
%             unaryImg = gray2rgb(255 * unaryPotential{ index });
%             segmentationImg = gray2rgb(255 * segmentation{ index });
            
            locationImg = getHeatmap(locationProbability{ index }, false);
            appearanceImg = getHeatmap(appearanceProbability{ index }, false);
            unaryImg = getHeatmap(unaryPotential{ index }, false);
            segmentationImg = gray2rgb(255 * segmentation{ index });

            fullframe = [frame, segImg; locationImg, appearanceImg; unaryImg, segmentationImg];
            writeVideo(writerObject, fullframe);
        end
        
    elseif(strcmp(mode, 'ShowProcessExperimental'))
        locationProbability = data.locationProbability;
        objectnessProbability = data.annotationProbability;
        appearanceProbability = data.appearanceProbability;
        unaryPotential = data.unaryPotential;
        for(index = 1: length(segmentation))
            frameid = options.ranges(range) + index - 1;
            frame = readFrame(options, frameid);

            segImg = overlaySegmentationMask(frame, ...
                segmentation{ index });

            locationImg = getHeatmap(locationProbability{ index }, false);
            objectnessImg = getHeatmap(objectnessProbability{ index }, false);
            appearanceImg = getHeatmap(appearanceProbability{ index }, false);
            unaryImg = getHeatmap(unaryPotential{ index }, false);
            
            fullframe = [ frame, segImg; ...
                locationImg, appearanceImg; ...
                objectnessImg, unaryImg ];
            writeVideo(writerObject, fullframe);
        end
    else
        error('Mode: "%s" not recognised', mode);
    end

    close(writerObject);
end

function result = overlaySegmentationBoundary(img, segmentation)

    se = strel('disk', 6, 0);
    [ height, width ] = size(segmentation);
    
    bounds = bwperim(segmentation);
    bounds = imdilate(bounds, se);
    bounds = bounds & (~segmentation);
    boundsindex = reshape(bounds, [], 1);
    linimg = reshape(img, [], 3);
    linimg(boundsindex, 1) = 0;
    linimg(boundsindex, 3) = 0;
    linimg(boundsindex, 2) = 255;
    
    result = reshape(linimg, height, width, 3);

end

function result = overlaySegmentationMask(img, segmentation)

    % Black background
    mask = gray2rgb(uint8(segmentation));
    result = img .* mask;
    
    % White background
%     [ height, width ] = size(segmentation);
%     indexes = reshape(~segmentation, [], 1);
%     result = reshape(img, [], 3);
%     result(indexes, :) = 255;
%     result = reshape(result, height, width, 3);

end

function result = overlayInMap(img, inMap)

    [ height, width ] = size(inMap);
        
    inmapIndex = reshape(inMap, [], 1);
    linimg = reshape(img, [], 3);
    linimg(inmapIndex, 1) = 255;
    linimg(inmapIndex, 2) = 0.25 * linimg(inmapIndex, 3);
    linimg(inmapIndex, 3) = 255;
    
    result = reshape(linimg, height, width, 3);

end

function result = overlayMotionBoundary(img, mbound)

    [ height, width ] = size(mbound);
        
    boundIndex = reshape(mbound, [], 1);
    linimg = reshape(img, [], 3);
    linimg(boundIndex, 1) = 255;
    linimg(boundIndex, 2) = 0;
    linimg(boundIndex, 3) = 0;
    
    result = reshape(linimg, height, width, 3);

end
