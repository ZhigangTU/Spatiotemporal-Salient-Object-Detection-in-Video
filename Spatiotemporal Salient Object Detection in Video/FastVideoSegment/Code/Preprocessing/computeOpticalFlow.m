% Wrapper to compute some given optical flow method
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

function flow = computeOpticalFlow( options, shot )
    
    flowfolder = fullfile( options.outfolder, 'flow', options.flowmethod );
    if( ~exist( flowfolder, 'dir' ) )
        mkdir( flowfolder );
    end

    fprintf( 'computeOpticalFlow: Processing shot #%i\n', shot );
    filename = fullfile( flowfolder, sprintf( 'flowShot%i.mat', shot ) );
    
    if( exist( filename, 'file' ) )
        % Shot already processed, skip
        fprintf( 'computeOpticalFlow: Data processed, skipping...\n' );
        return;
    else
        if( strcmp( options.flowmethod, 'broxPAMI2011' ) )
            flow = computeBroxPAMI2011Flow( options, shot );
        elseif( strcmp( options.flowmethod, 'sundaramECCV2010' ) )
            flow = computeSundaramECCV2010Flow( options, shot );
        else
            error( 'Unknown flow computation method' )
        end
        
        save( filename, 'flow', '-v7.3' );
    end
    
    fprintf( 'computeOpticalFlow: Shot #%i finished processing\n', shot );
    
end
