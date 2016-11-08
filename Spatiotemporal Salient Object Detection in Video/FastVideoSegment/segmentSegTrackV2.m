% Script that runs the segmentation algorithm on the SegTrack dataset
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

clear
addpath( genpath( '.' ) )

foldername = fileparts( mfilename( 'fullpath' ) );

options.infolder = fullfile( foldername, 'Data', 'inputs', 'SegTrackV2' );
options.outfolder = fullfile( foldername, 'Data', 'outputs', 'SegTrackV2' );
options.flowmethod = 'broxPAMI2011';
options.superpixels = 'Turbopixels';
options.maxedge = inf;
options.vocal = true;
options.visualise = false;
options.ranges = [ 1, 99, 141, 177, 251, 530, 559, 590, 622, 865 ];

videoid = { 'bird_of_paradise', 'penguin', 'bmx', 'drift', 'frog',...
    'hummingbird', 'monkey', 'soldier', 'worm' };

segmfolder = fullfile( options.outfolder, 'segmentations', ...
    'VideoRapidSegment' );
if( ~exist( segmfolder, 'dir' ) ), mkdir( segmfolder ), end;

for( shot = 9: length( options.ranges ) - 1 )
   
    data.flow = loadFlow( options, shot );
    if( isempty( data.flow ) )
        data.flow = computeOpticalFlow( options, shot );
    end
    
    data.superpixels = loadSuperpixels( options, shot );
    if( isempty( data.superpixels ) )
        data.superpixels = computeSuperpixels( options, shot );
    end
    
    data.imgs = readAllFrames( options, shot );
    data.id = shot;
    
    params = getSegTrackParamsV2( videoid{ shot } );
    
    
    if ~exist([options.infolder '\GroundTruth\groundTruthShot' num2str(shot) '.mat'],'file')
        
        startFrame = options.ranges(shot);
        endFrame = options.ranges(shot+1) - 1;
        groundTruth = cell((endFrame-startFrame+1),1);
        for grni = 1 : (endFrame - startFrame + 1)
            groundTruth{grni} = imread([options.infolder '\GroundTruth\' sprintf( ...
                '%08d.png', grni + startFrame - 1 )]); 
            if size(groundTruth{grni},3) == 3 
                groundTruth{grni} = rgb2gray(groundTruth{grni})~=0;
            end
        end


        save([options.infolder '\GroundTruth\groundTruthShot' num2str(shot) '.mat'],...
            'groundTruth', '-v7.3');
    end
    
    load([options.infolder '\GroundTruth\groundTruthShot' num2str(shot) '.mat'])
    segmentation = videoRapidSegment( options, params, data );
    segmentation = getLargestSegmentAndNeighbours( segmentation , params );
  
    
    filename = fullfile( segmfolder, ...
        sprintf( 'segmentation-%s.mat', videoid{ shot } ) );
    save( filename, 'segmentation', '-v7.3' );
    
    avgMislabelled = ...
        getAverageMislabelledPixels_V2( options, shot, segmentation );
    
    figure;
    subplot(2,1,1); imshow(segmentation{1});
    subplot(2,1,2); imshow(groundTruth{1});
    supertitle(['Loc-',num2str(params.locationWeight),' Temp-'...
    ,num2str(params.temporalWeight),' Spa-',num2str(params.spatialWeight)...
    ,' Average Error:', num2str(avgMislabelled)])
    
    
    fprintf( 'Average number of mislabelled pixels for %s: %i\n', ...
        videoid{ shot }, avgMislabelled );
    
end

rmpath( genpath( '.' ) )
