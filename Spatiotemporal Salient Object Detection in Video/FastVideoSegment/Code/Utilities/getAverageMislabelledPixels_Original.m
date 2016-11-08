% Funtion that counts the average number of mislabelled pixels according to
% the ground truth
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

function avgMislabelled = getAverageMislabelledPixels( options, shot, segmentation )
    
    gtFile = fullfile( options.infolder, 'GroundTruth', sprintf( 'groundTruthShot%i.mat', shot ) );
    groundTruth = load( gtFile );
    groundTruth = groundTruth.groundTruth;

    frames = length( groundTruth );
    
    totalTrue = 0;
    totalPredicted = 0;
    totalCorrect = 0;
    for frame = 1: frames 

        [ truePoints, predictedPoints, correctPoints ] = ...
            compareToGroundTruth( segmentation{ frame }, ...
            groundTruth{ frame } );
            
        totalTrue = totalTrue + truePoints;
        totalPredicted = totalPredicted + predictedPoints;
        totalCorrect = totalCorrect + correctPoints;

    end
    
    wrong = totalTrue - totalCorrect + totalPredicted - totalCorrect;
    avgMislabelled = round( wrong / frames );
    
end
