% Function to create a grayscale image based on given 2-label potentials
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

function result = mapBelief( potential )

    if( iscell( potential ) )
        frames = length( potential );
        result = cell( frames, 1 );
        for( i = 1: frames )
            result{ i } = mapBelief( potential{ i } );
        end
    else
        [ ~, ~, ndims ] = size( potential );
        if( ndims == 1 )
            potentialSum = sum( potential, 2 );
            result = 1 - potential( :, 1 ) ./ potentialSum;
        else
            potentialSum = potential( :, :, 1 ) + potential( :, :, 2 );
            result = 1 - potential( :, :, 1 ) ./ potentialSum;
            
            infindsNom = isinf( potential( :, :, 1 ) );
            infindsDen = isinf( potential( :, :, 2 ) );
            
            result( infindsNom ) = 0;
            result( infindsNom & infindsDen ) = 0.5;
            
        end
        
    end

end
