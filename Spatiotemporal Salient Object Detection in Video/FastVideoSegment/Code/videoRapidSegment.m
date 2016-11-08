function segmentation = videoRapidSegment(options, params, data)
    
    % Params parsing
    if(~isfield(params, 'locationWeight') || isempty(params.locationWeight))
        params.locationWeight = 1;
    end
    
    if(~isfield(params, 'spatialWeight') || isempty(params.spatialWeight))
        params.spatialWeight = 5000;
    end
    
    if(~isfield(params, 'temporalWeight') || isempty(params.temporalWeight))
        params.temporalWeight = 4000;
    end
    
    if(~isfield(params, 'fadeout') || isempty( params.fadeout))
        params.fadeout = 0.0001;
    end
    
    if(~isfield(params, 'maxIterations') || isempty( params.maxIterations))
        params.maxIterations = 4;
    end
    
    if( isfield(params, 'foregroundMixtures') && ~isempty( params.foregroundMixtures))
        fgMix = params.foregroundMixtures;
    else
        fgMix = 5;
    end
    
    if(isfield(params, 'backgroundMixtures' ) && ~isempty( params.backgroundMixtures))
        bgMix = params.backgroundMixtures;
    else
        bgMix = 8;
    end
    
    if(isfield(params, 'locationNorm' ) && ~isempty( params.locationNorm))
        locationNorm = params.locationNorm;
    else
        locationNorm = 0.75;
    end
    
    if(~isfield( options, 'visualise') || isempty(options.visualise))
        options.visualise = false;
    end
    if(~isfield( options, 'vocal') || isempty(options.vocal))
        options.vocal = false;
    end
    % End of params parsing
    
    Images  = data.imgs;
    magflow = data.flowmag;
    Orgflow = data.Orgflow;
    flow    = data.flow;   %flow--int16
    superpixels = data.superpixels;
    
    % Compute inside-outside maps
    if( options.vocal )
        tic; 
        fprintf( 'videoRapidSegment: Computing inside-outside maps...\t' ); 
    end
    
    data.inMaps = getInOutMaps(Orgflow, flow, magflow, Images);    % Orgflow; flow
    inRatios = getSuperpixelInRatio(superpixels, data.inMaps);
    if options.vocal
        toc; 
    end
    
    imgs = data.imgs;
    frames = length(flow);
    
    % Compute pairwise potentials
    if options.vocal
        tic; 
        fprintf( 'videoRapidSegment: Computing pairwise potentials...\t' ); 
    end
    
    [superpixels, nodeFrameId, bounds, labels] = makeSuperpixelIndexUnique(superpixels);
    [colours, centres, ~] = getSuperpixelStats(imgs, superpixels, labels);
    
    pairPotentials = computePairwisePotentials(params, superpixels, flow, colours, centres, labels);
    
    if options.vocal
        toc; 
    end
    
    colours = uint8(round(colours));
    
    % Preallocate space for unary potentials
    nodes = size( colours, 1 );
    unaryPotentials = zeros( nodes, 2 );

    % Create location priors
    if options.vocal
        tic; 
        fprintf( 'videoRapidSegment: Computing location priors...\t\t' ); 
    end
    [~, accumulatedInRatios] = accumulateInOutMap(params, data);
    locationMasks = cell2mat(accumulatedInRatios);
    
    locationUnaries = 0.5 * ones( nodes, 2, 'single' );

    locationUnaries(1: length(locationMasks), 1) = locationMasks / (locationNorm * max(locationMasks));
    locationUnaries(locationUnaries > 0.95) = 0.999;
    
    for frame = 1: frames
        start = bounds( frame );
        stop = bounds( frame + 1 ) - 1;
        
        frameMasks = locationUnaries(start: stop, 1);
        overThres = sum(frameMasks > 0.6) / single((stop - start));

        if overThres < 0.05
            E = 0.005;
        else
            E = 0.000;
        end
        locationUnaries(start: stop, 1) = max(locationUnaries(start: stop, 1), E);      
    end
    locationUnaries(:, 2) = 1 - locationUnaries(:, 1);
    
    if options.vocal
        toc; 
    end

    masks = 0.19 * ones( nodes, 1 );
    masks( 1: bounds( frames + 1 ) - 1 ) = single( cell2mat( inRatios ) );

    % Create binary masks for foreground/background initialisation
    foregroundMasks = masks > 0.2;
    backgroundMasks = masks < 0.05;

    % Create fading frame weight
    weights = zeros( 1 + 2 * frames, 1, 'single' );
    middle = frames + 1;
    for i = 1: length( weights ) 
        weights( i ) = exp( - params.fadeout * ( i - middle ) ^ 2 );
    end

    % Initialise segmentations
    if options.vocal
        tic; fprintf( 'videoRapidSegment: Computing initial segmentation...\t' ); 
    end
    fgColors = colours( foregroundMasks, : );
    bgColors = colours( backgroundMasks, : );
    for frame = 1: frames        
        ids = nodeFrameId - frame + middle;
        
        fgNodeWeights = masks(foregroundMasks).*weights(ids(foregroundMasks));
        bgNodeWeights = (1 - masks(backgroundMasks)).*weights(ids(backgroundMasks));

        [uniqueFgColours, fgNodeWeights] = findUniqueColourWeights(fgColors, fgNodeWeights);
        [uniqueBgColours, bgNodeWeights] = findUniqueColourWeights(bgColors, bgNodeWeights);
        
        startIndex = bounds(frame);
        stopIndex  = bounds(frame + 1) - 1;
        
        if( size( uniqueFgColours, 1 ) < fgMix || size( uniqueBgColours, 1 ) < bgMix )
            warning( 'Too few data points to fit GMM...\n' ); %#ok<WNTAG>
            unaryPotentials( startIndex: stopIndex, : ) = -log( 0.5 );
        else
            [ fgModel ] = fitGMM( fgMix, uniqueFgColours, fgNodeWeights );
            [ bgModel ] = fitGMM( bgMix, uniqueBgColours, bgNodeWeights );

            appearanceUnary = getUnaryAppearance(single(colours(nodeFrameId == frame, :)), fgModel, bgModel);

            unaryPotentials(startIndex: stopIndex, :) = -params.locationWeight * log(locationUnaries(startIndex: stopIndex, :)) + -log(appearanceUnary);
        end
    end

    [~, labels] = maxflow_mex_optimisedWrapper(pairPotentials, single(10 * unaryPotentials));
    
    segmentation = superpixelToPixel(labels, superpixels);
    
    if(options.vocal), toc; end
    
    % Check that we did not get a trivial, all-background/foreground segmentation
    if( all( labels ) || all( ~labels ) )
        if( options.vocal )
            fprintf( 'videoRapidSegment: Trivial segmentation detected, exiting...\n' ); 
        end
        return;
    end
    
    % Iterating segmentations
    for i = 2: params.maxIterations 
        if( options.vocal ), tic; fprintf( 'videoRapidSegment: Iteration: %i...\t\t\t', i ); end
        
        fgColors = colours( labels, : );
        bgColors = colours( ~labels, : );
        for frame = 1: frames
            oldLabels = labels;

            ids = nodeFrameId - frame + middle;

            fgNodeWeights = weights( ids( labels ) );
            bgNodeWeights = weights( ids( ~labels ) );

            [ uniqueFgColours, fgNodeWeights ] = ...
                findUniqueColourWeights( fgColors, fgNodeWeights );
            [ uniqueBgColours, bgNodeWeights ] = ...
                findUniqueColourWeights( bgColors, bgNodeWeights );

            if( size( uniqueFgColours, 1 ) < fgMix || ...
                size( uniqueBgColours, 1 ) < bgMix )
                warning( 'videoRapidSegment: Too few data points to fit GMM...\n' ); %#ok<WNTAG>
                return;
            end
            
            [ fgModel ] = fitGMM( fgMix, uniqueFgColours, fgNodeWeights );
            [ bgModel ] = fitGMM( bgMix, uniqueBgColours, bgNodeWeights );
            
            appearanceUnary = getUnaryAppearance( ...
                single( colours( nodeFrameId == frame, : ) ), ...
                fgModel, bgModel );

            startIndex = bounds( frame );
            stopIndex = bounds( frame + 1 ) - 1;

            unaryPotentials( startIndex: stopIndex, : ) = ...
                -params.locationWeight * log( ...
                locationUnaries( startIndex: stopIndex, : ) ) + ...
                -log( appearanceUnary );
        end
        
        [ ~, labels ] = maxflow_mex_optimisedWrapper( pairPotentials, single( 10 * unaryPotentials ) );
        segmentation = superpixelToPixel( labels, superpixels );
        
        if( options.vocal ), toc; end
        
        if( ( i == params.maxIterations ) || all( all( oldLabels == labels ) ) )

            if( options.vocal )
                fprintf( 'videoRapidSegment: Convergence or maximum number of iterations reached\n' ); 
            end

            if( options.visualise )
                if( options.vocal )
                    tic; 
                    fprintf( 'videoRapidSegment: Creating segmentation video...\t' ); 
                end

                videoParams.name = sprintf( 'segmentation%i', data.id );
                videoParams.range = data.id;
                mode = 'ShowProcess';

                data.locationProbability = superpixelToPixel( ...
                    locationUnaries( :, 1 ), superpixels );
                data.appearanceProbability = mapBelief( ...
                    superpixelToPixel( - log( single( ...
                    getUnaryAppearance( single( colours ), fgModel, ...
                    bgModel ) ) ), superpixels ) );
                data.unaryPotential = mapBelief( superpixelToPixel( ...
                    unaryPotentials, superpixels ) );

                data.segmentation = segmentation;

                createSegmentationVideo( options, videoParams, data, mode );

                videoParams.name = sprintf( 'segmentation%i-dominantObject', data.id );
                data.segmentation = getLargestSegmentAndNeighbours(segmentation);
                createSegmentationVideo( options, videoParams, data, mode );
                clear data

                if( options.vocal ), toc; end
            end
            
            break;
        end
        
        % Check that we did not get a trivial, all-background/foreground segmentation
        if( all( labels ) || all( ~labels ) )
            if options.vocal
                fprintf( 'videoRapidSegment: Trivial segmentation detected, exiting...\n' ); 
            end
            return;
        end
        
    end

    if options.vocal
        fprintf( 'videoRapidSegment: Algorithm stopped after %i iterations.\n', i ); 
    end

end

function potentials = computePairwisePotentials( params, superpixels, flow, colours, centres, labels )   % flow--int16

    [sSource, sDestination] = getSpatialConnections( superpixels, labels );
    [tSource, tDestination, tConnections] = getTemporalConnections(flow, superpixels, labels);

    sSqrColourDistance = sum( ( colours( sSource + 1, : ) - colours( sDestination + 1, : ) ) .^ 2, 2 ) ;
    sCentreDistance = sqrt( sum( ( centres( sSource + 1, : ) - centres( sDestination + 1, : ) ) .^ 2, 2 ) );
    tSqrColourDistance = sum( ( colours( tSource + 1, : ) - colours( tDestination + 1, : ) ) .^ 2, 2 ) ;

    sBeta = 0.5 / mean( sSqrColourDistance ./ sCentreDistance );
    tBeta = 0.5 / mean( tSqrColourDistance .* tConnections );

    sWeights = exp(-sBeta * sSqrColourDistance) ./ sCentreDistance;
    tWeights = tConnections .* exp(-tBeta * tSqrColourDistance);

    potentials.source = [sSource; tSource];
    potentials.destination = [sDestination; tDestination];
    potentials.value = [params.spatialWeight * sWeights; params.temporalWeight * tWeights];
    
end
