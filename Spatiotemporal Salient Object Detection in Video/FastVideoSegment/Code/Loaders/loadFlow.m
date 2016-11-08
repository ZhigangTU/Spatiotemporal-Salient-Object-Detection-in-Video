function result = loadFlow(flowDir, range, flowmethod)

    file = fullfile(flowDir, flowmethod, sprintf('flowShot%i.mat', range));
    if( exist( file, 'file' ) )
        flow = load( file );
        if( isfield( flow, 'input' ) )
            result = flow.input;
        elseif( isfield( flow, 'flow' ) )
            result = flow.flow;
        else
            warning( '%s: no known field found\n', file );
            result = [];
        end
    else
        warning( '%s not found\n', file );
        result = [];
    end

end
