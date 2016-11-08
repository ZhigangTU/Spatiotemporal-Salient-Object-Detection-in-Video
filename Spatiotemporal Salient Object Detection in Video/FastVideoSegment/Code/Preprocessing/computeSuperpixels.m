% Function to compute some given superpixel method for a given shot
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

function superpixels = computeSuperpixels( options, shot, numImages, framesdata)
    
    superpixelfolder = fullfile( options.SPoutfolder, 'superpixels', options.superpixels );
    if( ~exist( superpixelfolder, 'dir' ) )
        mkdir( superpixelfolder );
    end

    fprintf( 'computeSuperpixels: Processing shot #%i\n', shot );
    filename = fullfile(superpixelfolder, sprintf( 'superpixelsShot%i.mat', shot ) );
    
    if( exist( filename, 'file' ) )
        % Shot already processed, skip
        fprintf( 'computeSuperpixels: Data processed, skipping...\n' );
        superpixels = loadSuperpixels(options, shot);
        return;
    else
        if( strcmp( options.superpixels, 'Turbopixels' ) )
            superpixels = computeTurbopixels( options, shot );
        elseif( strcmp( options.superpixels, 'SLIC' ) )
            superpixels = computeSLIC(options, numImages, framesdata);
        else
            error( 'Unknown superpixel oversegmentation method' )
        end
        
        save( filename, 'superpixels', '-v7.3' );
    end
    
    fprintf( 'computeSuperpixels: Shot #%i finished processing\n', shot );
    
end
