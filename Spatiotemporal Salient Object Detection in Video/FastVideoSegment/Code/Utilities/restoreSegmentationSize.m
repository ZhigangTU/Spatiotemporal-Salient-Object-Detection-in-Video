function segmentation = restoreSegmentationSize( options, shot, ...
    segmentation )

    options.maxedge = inf;
    [ oHeight, oWidth, ~ ] = ...
        size( readFrame( options, options.ranges( shot ) ) );  
    
    if( iscell( segmentation ) )
        [ height, width ] = size( segmentation{ 1 } );
        if( height ~= oHeight || width ~= oWidth )
            frames = length( segmentation );
            for( frame = 1: frames )
                segmentation{ frame } = ...
                    imresize( segmentation{ frame }, [ oHeight, oWidth ] );
            end
        end
    else
        [ height, width ] = size( segmentation );
        if( height ~= oHeight || width ~= oWidth )
            segmentation = imresize( segmentation, [ oHeight, oWidth ] );
        end
    end

end
