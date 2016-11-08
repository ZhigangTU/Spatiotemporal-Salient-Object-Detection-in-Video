% Function to produce inside-outside maps of a shot given the optical flow
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

function output = getInOutMaps( flow )

    frames = length( flow );
    output = cell( frames, 1 );

    [ height, width, ~ ] = size( flow{ 1 } );
    
    % Motion boundaries touching the edges will be cut!
    sideCut = false( height, width );
    sideCut( 1: 20, : ) = true;
    sideCut( end - 20: end, : ) = true;
    sideCut( :, 1: 20 ) = true;
    sideCut( :, end - 20: end ) = true;
    
    for( frame = 1: frames )
        boundaryMap = getProbabilityEdge( flow{ frame }, 3 );

        inVotes = getInPoints( boundaryMap, sideCut, false );
        
        if( getFrameQuality( inVotes > 4 ) < 0.2 )
            boundaryMap = calibrateProbabilityEdge( flow{ frame }, 0.71 );
            inVotes = getInPoints( boundaryMap, sideCut, false );
        end
        
        output{ frame } = inVotes > 4;
    end    
    
end
