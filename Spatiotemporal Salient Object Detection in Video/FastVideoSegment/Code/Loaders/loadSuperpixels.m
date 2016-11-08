function superpixels = loadSuperpixels(options, flowPath, shotid)
    file = fullfile(flowPath, 'superpixels', options.superpixels, sprintf( 'superpixelsShot%i.mat', shotid));
    if( exist( file, 'file' ) )
        superpixels = load( file );
        if( isfield( superpixels, 'input' ) )
            superpixels = superpixels.input;
        elseif( isfield( superpixels, 'superpixels' ) )
            superpixels = superpixels.superpixels;
        elseif( isfield( superpixels, 'superPixels' ) )
            superpixels = superpixels.superPixels;
        else
        		warning( '%s: no known field found\n', file );
            superpixels = [];
        end
    else
        warning( '%s not found\n', file );
        superpixels = [];
    end
    
end