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

addpath( genpath( '.' ) )

foldername = fileparts( mfilename( 'fullpath' ) );

options.infolder = fullfile( foldername, 'Data', 'inputs', 'SegTrack' );
options.outfolder = fullfile( foldername, 'Data', 'outputs', 'SegTrack' );
options.flowmethod = 'broxPAMI2011';
options.superpixels = 'Turbopixels';
options.maxedge = inf;
options.vocal = false;
options.visualise = true;
options.ranges = [ 1, 31, 60, 81, 152, 203 ];

videoid = { 'birdfall', 'cheetah', 'girl', 'monkey', 'parachute' };

segmfolder = fullfile( options.outfolder, 'segmentations', ...
    'VideoRapidSegment' );
if( ~exist( segmfolder, 'dir' ) ), mkdir( segmfolder ), end;

for( shot = 1: length( options.ranges ) - 1 )
   
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
    
    params = getSegTrackParams( videoid{ shot } );
        
    segmentation = videoRapidSegment( options, params, data );
    segmentation = getLargestSegmentAndNeighbours( segmentation );
    filename = fullfile( segmfolder, ...
        sprintf( 'segmentation-%s.mat', videoid{ shot } ) );
    save( filename, 'segmentation', '-v7.3' );
    
    avgMislabelled = ...
        getAverageMislabelledPixels( options, shot, segmentation );
    fprintf( 'Average number of mislabelled pixels for %s: %i\n', ...
        videoid{ shot }, avgMislabelled );
    
end

rmpath( genpath( '.' ) )
