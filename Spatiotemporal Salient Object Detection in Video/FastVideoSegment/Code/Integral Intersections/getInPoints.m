% Returns a matrix representing the number of rays that "believe" that each
% pixel is inside an object
%
% Inputs:
%   edgeMap: A 2-D matrix representing the belief that a pixel is a
%       boundary.
%   cutMask: A boolean mask of pixels. Boundary lines overlapping with the
%       mask are removed (belief of pixel being a boundary is set to 0).
%       Default cutMask: Pixels that are less than 21 pixels away from the
%       matrix borders are set to true.
%   thresholded: Boolean variable, indicating whether the edgeMap is
%       already thresholded. If true, any pixel in the edgeMap with belief
%       greater than 0 is considered to be an edge. If false, any pixel
%       with belief greater than 0.5 is considered to be an edge.
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

function result = getInPoints( edgeMap, cutMask, thresholded )

    if( ~exist( 'cutMask', 'var' ) || isempty( cutMask ) )
        cutMask = false( size( edgeMap ) );
        cutMask( 1: 20, : ) = true;
        cutMask( end - 20: end, : ) = true;
        cutMask( :, 1: 20 ) = true;
        cutMask( :, end - 20: end ) = true;
    end
    
    if( ~exist( 'thresholded', 'var' ) || isempty( thresholded ) || ...
        ~thresholded )
    
        mask = zeros( size( edgeMap ), 'uint8' );
        mask( edgeMap > 0.5 ) = 1;
    else
        mask = edgeMap;
    end

    mask = uint8( cleanBoundaries( logical( mask ), cutMask ) );

    result = integralIntersections( mask );
    % A point on a boundary cannot be "in"
    result( logical( mask ) ) = 0;

end
