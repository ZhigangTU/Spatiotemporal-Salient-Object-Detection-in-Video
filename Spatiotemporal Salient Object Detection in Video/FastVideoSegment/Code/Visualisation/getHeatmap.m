% Function to create a heat map from a greyscale map.
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

function heatmap = getHeatmap( map, normalise )

    colors = uint8( [ 0, 0, 143; 0, 0, 159; 0, 0, 175; 0, 0, 191; 0, 0, 207; ...
        0, 0, 223; 0, 0, 239; 0, 0, 255; 0, 16, 255; 0, 32, 255; ...
        0, 48, 255; 0, 64, 255; 0, 80, 255; 0, 96, 255; 0, 112, 255; ...
        0, 128, 255; 0, 143, 255; 0, 159, 255; 0, 175, 255; 0, 191, 255; ...
        0, 207, 255; 0, 223, 255; 0, 239, 255; 0, 255, 255; ...
        16, 255, 239; 32, 255, 223; 48, 255, 207; 64, 255, 191; ...
        80, 255, 175; 96, 255, 159; 112, 255, 143; 128, 255, 128; ...
        143, 255, 112; 159, 255, 96; 175, 255, 80; 191, 255, 64; ...
        207, 255, 48; 223, 255, 32; 239, 255, 16; 255, 255, 0; ...
        255, 239, 0; 255, 223, 0; 255, 207, 0; 255, 191, 0; ...
        255, 175, 0; 255, 159, 0; 255, 143, 0; 255, 128, 0; ...
        255, 112, 0; 255, 96, 0; 255, 80, 0; 255, 64, 0; 255, 48, 0; ...
        255, 32, 0; 255, 16, 0; 255, 0, 0; 239, 0, 0; 223, 0, 0; ...
        207, 0, 0; 191, 0, 0; 175, 0, 0; 159, 0, 0; 143, 0, 0; 128, 0, 0 ] );
    
    if( ~exist( 'normalise', 'var' ) || isempty( normalise ) )
    	normalise = true;
  	end
    
    [ height, width ] = size( map );
    
    if( normalise )
    
	  minmap = minmin( map );
  	  maxmap = maxmax( map );
      %disp([minmap, maxmap])
  	  colormapping = round( 63 * ( map - minmap ) / ( maxmap - minmap ) ) + 1;

    else
      colormapping = round( 63 * map ) + 1;
      %disp([minmin(colormapping), maxmax(colormapping)])
    end

    % Sanity checks
    colormapping( isnan( colormapping ) ) = 32;
    colormapping( isinf( colormapping ) ) = 32;

    heatmap = colors( reshape( colormapping, [], 1 ), : );
    heatmap = reshape( heatmap, height, width, 3 );

end
