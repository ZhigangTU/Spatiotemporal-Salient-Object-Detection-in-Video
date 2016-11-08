function avgMislabelled = getAverageMislabelledPixels( options, segmentation )   
%     gtFile = fullfile( options.infolder, 'GroundTruth', sprintf( 'groundTruthShot%i.mat', shot ) );
%     groundTruth = load( gtFile );
%     groundTruth = groundTruth.groundTruth;
%     gtFile = fullfile(options.infolder, 'GroundTruth'); %'1'
%     groundTruth = readAllFrames(gtFile);
    
    gtFile = fullfile(options.infolder, 'GroundTruth', sprintf('GT%s.mat', options.SeqName)); %'1'
    groundTruth = load( gtFile );
    groundTruth = groundTruth.groundTruth;
    
    frames = length(groundTruth);    
    totalTrue = 0;
    totalPredicted = 0;
    totalCorrect = 0;
    for frame = 1: frames 
        [truePoints, predictedPoints, correctPoints] = compareToGroundTruth(segmentation{frame}, groundTruth{frame});            
        totalTrue = totalTrue + truePoints;
        totalPredicted = totalPredicted + predictedPoints;
        totalCorrect = totalCorrect + correctPoints;
    end
    
    wrong = totalTrue - totalCorrect + totalPredicted - totalCorrect;
    avgMislabelled = round( wrong / frames );
    
end
