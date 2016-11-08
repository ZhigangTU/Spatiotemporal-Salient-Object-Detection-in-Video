% Function to find the largest connected component in each frame, and very close
% neighbouring connected components
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

function output = getLargestSegmentAndNeighbours( map, params )

    % params.ratio contains the ratio of the radius of the larsgest
    % segments equivalent area circle by which the largest segment gets
    % expanded
    if( exist( 'params', 'var' ) && ...
        isfield( params, 'ratio' ) && ...
        ~isempty( params.ratio ) )
        ratio = params.ratio;
    else
        ratio = 0.3;
    end

    if( ~iscell( map ) )
        output = getLargetSegmentAndNeighboursInFrame( map, ratio );
    else
        output = cell( size( map ) );
        
        for( i = 1: length( map ) )
            output{ i } = getLargetSegmentAndNeighboursInFrame( ...
                map{ i }, ratio );
        end
    end
    
end

function output = getLargetSegmentAndNeighboursInFrame( map, ratio )

    % Get connected components
    [ labelmap, labels ] = bwlabel( map, 8 );
    
    % If there are no components, just return
    if( labels == 0 )
        output = map;
        return; 
    end
    
    % Get largest segment
    segments = regionprops( labelmap, 'area' );
    largest = 0;
    largestSize = 0;
    for( i = 1: labels )
        if( segments(i).Area > largestSize )
            largestSize = segments(i).Area;
            largest = i;
        end
    end
    
    % Expand largest segment
    best = labelmap == largest;
    disk = strel( 'disk', round( ratio * sqrt( largestSize / pi ) ), 0);
    best = imdilate( best, disk );
    
    % Connect expanded segment to ovelapping segments
    best = best | map;
    
    % Find new largest segment
    [ expanded, labels ] = bwlabel( best, 8 );
    segments = regionprops( expanded, 'area' );
    largest = 0;
    largestSize = 0;
    for( i = 1: labels )
        if( segments(i).Area > largestSize )
            largestSize = segments(i).Area;
            largest = i;
        end
    end
    
    % Return the intersection between largest expanded segment and original
    output = ( expanded == largest ) & map;

end
