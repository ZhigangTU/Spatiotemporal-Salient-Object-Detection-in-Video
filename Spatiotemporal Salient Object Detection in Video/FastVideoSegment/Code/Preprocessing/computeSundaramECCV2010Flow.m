% Wrapper to compute the sundaramECCV2010 optical flow method for a given shot
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

function flow = computeSundaramECCV2010Flow( options, shot )

    startFrame = options.ranges( shot );
    endFrame = options.ranges( shot + 1 ) - 1;
    flowframes = endFrame - startFrame;
    
    flow = cell( 1, flowframes );

    totalTimeTaken = 0;
    
    currImage = readFrame( options, startFrame );
    for( i = startFrame + 1: endFrame )
        tic

        index = i - startFrame;
        nextImage = readFrame( options, i );

        if( options.vocal )
            fprintf( 'Computing optical flow of pair: %i of %i... ', ...
                index, flowframes );
        end

        [ forwardFlow, ~ ] = ...
            sundaramECCV10_ldof_GPU_mex( currImage, nextImage );
        flow{ index } = int16( forwardFlow );
        
        currImage = nextImage;
        timeTaken = toc;
        totalTimeTaken = totalTimeTaken + timeTaken;
        
        if( options.vocal )
            fprintf( 'done. Time taken: %.2f sec\n', timeTaken );
        end

    end
    
    if( options.vocal )
        fprintf( 'Total time taken: %.2f sec\n', totalTimeTaken );
        fprintf( 'Average time taken per frame: %.2f sec\n', ...
            totalTimeTaken / flowframes );
    end
    
end
