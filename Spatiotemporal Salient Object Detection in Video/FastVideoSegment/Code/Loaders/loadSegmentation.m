% Function to load the computed VRS segmentation
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

function result = loadSegmentation( options, range )

    file = fullfile( options.outfolder, 'segmentations', ...
        'VideoRapidSegment', sprintf( 'segmentationShot%i.mat', range ) );
    if( exist( file, 'file' ) )
        input = load( file );
        if( isfield( input, 'input' ) )
            result = input.input;
        elseif( isfield( input, 'segmentation' ) )
            result = input.segmentation;
        else
            result = [];
        end
    else
        result = [];
    end

end
