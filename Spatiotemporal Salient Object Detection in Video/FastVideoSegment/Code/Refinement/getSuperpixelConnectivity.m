% Function to compute the connections between superpixels given the flow
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

function connectivity = getSuperpixelConnectivity( flow, superPixels )

    frames = length( flow );
    connectivity = cell( frames, 1 );
    for( frame = 1: frames )
        sPixels = double( maxmax( superPixels{ frame } ) );
        sPixelsNext = double( maxmax( superPixels{ frame + 1 } ) );

        connectivity{ frame } = sparse( double( ...
            superPixelConnectivity( flow{ frame }, ...
            superPixels{ frame }, superPixels{ frame + 1 }, ...
            sPixels, sPixelsNext ) ) );
    end

end
