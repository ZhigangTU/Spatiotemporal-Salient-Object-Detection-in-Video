function [ truePoints, labeled, correctlyLabeled ] = ...
    compareToGroundTruth( input, groundTruth )

    truePoints = sum( sum( groundTruth ) );
    labeled = sum( sum( input ) );
    correctlyLabeled = sum( sum( input & groundTruth ) );

end