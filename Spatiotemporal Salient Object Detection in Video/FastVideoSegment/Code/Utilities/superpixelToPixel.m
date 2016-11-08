% Function to label pixels given their superpixel ids
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

function result = superpixelToPixel( labels, superpixels )

    frames = length( superpixels );
    ndims = size( labels, 2 );
    
    result = cell( frames, 1 );
    
    if( ndims == 1 )
        for( i = 1: frames )
            result{ i } = labels( double( superpixels{ i } ) );
        end
    else
        for( i = 1: frames )
            [ height, width ] = size( superpixels{ i } );
            spixels = double( reshape( superpixels{ i }, [], 1 ) );
            
            result{ i } = labels( spixels, : );
            result{ i } = reshape( result{ i }, height, width, ndims );
        end
    end

end
