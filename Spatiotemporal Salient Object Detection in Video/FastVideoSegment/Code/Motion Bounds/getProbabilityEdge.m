% Function to compute motion boundary probabilities base on optical flow
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

function result = getProbabilityEdge( flowframe, mode, gradLambda )

    if( ~exist( 'mode', 'var' ) || isempty( mode ) )
        mode = 3;
    end
    
    if( ischar( mode ) )
        if( strcmp( mode, 'gradient' ) )
            mode = 1;
        elseif( strcmp( mode, 'orientation' ) )
            mode = 2;
        elseif( strcmp( mode, 'gradient+orientation' ) )
            mode = 3;
        else
            error( 'Unknown or invalid mode selected' );
        end
    end

    if( ~exist( 'gradLambda', 'var' ) || isempty( gradLambda ) )
        gradLambda = 0.7;
    end
    
    if( mode == 1 )
        gradient = getFlowGradient( flowframe );
        magnitude = getMagnitude( gradient );

        result = 1 - exp( -gradLambda * magnitude );
    elseif( mode == 2 )
        result = getFlowDifference( flowframe );
    elseif( mode == 3 )
        gradient = getFlowGradient( flowframe );
        magnitude = getMagnitude( gradient );
        
        gradBoundary = 1 - exp( -gradLambda * magnitude );        
        rotBoundary = getFlowDifference( flowframe );
        
        large = gradBoundary > 0.6 ;
        medium = gradBoundary <= 0.6 & gradBoundary > 0.25;
        result = 0.1 * gradBoundary;
        result( large ) = gradBoundary( large );
        result( medium ) = ( gradBoundary( medium ) .* rotBoundary( medium ) );
    else
        error( 'Unknown or invalid mode selected' );
    end
    
end
