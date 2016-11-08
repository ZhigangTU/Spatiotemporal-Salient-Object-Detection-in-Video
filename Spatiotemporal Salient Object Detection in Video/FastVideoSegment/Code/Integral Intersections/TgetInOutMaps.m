function output = TgetInOutMaps(flow, magflow, Images, Orgflow)
    frames = length( flow );
    FInOut = cell( frames, 1 );
    OrgFInOut = cell( frames, 1 );

    [ height, width, ~ ] = size( flow{ 1 } );
    
    % Motion boundaries touching the edges will be cut!
    sideCut = false( height, width );
    sideCut( 1: 20, : ) = true;
    sideCut( end - 20: end, : ) = true;
    sideCut( :, 1: 20 ) = true;
    sideCut( :, end - 20: end ) = true;
    
%     % Method 1: Original Method
%     for frame = 1:frames
%         boundaryMap = getProbabilityEdge( flow{ frame }, 3 );  % getProbabilityEdge(flowframe, mode, gradLambda)
% 
%         inVotes = getInPoints( boundaryMap, sideCut, false );
%         
%         if( getFrameQuality( inVotes > 4 ) < 0.2 )
%             boundaryMap = calibrateProbabilityEdge( flow{ frame }, 0.71 );
%             inVotes = getInPoints( boundaryMap, sideCut, false );
%         end
%         
%         output{ frame } = inVotes > 4;
%     end
    
    % Method 2: Learning to Detect Motion Boundaries
    Tesh = 0.20;
    if isempty(Orgflow)
        IndxSmallFlowS1 = [];
        for frame = 1:frames
            IndxImg = uint8(Images{frame});
            Indxmagflow = magflow{frame};

            Indxflow    = flow{frame};

            boundary_ColorMap     = detect_motionboundaries(IndxImg);      % not accurate
            boundary_ColorFlowMap = detect_motionboundaries(IndxImg, Indxflow);

            if max(Indxmagflow(:)) > 1.0
                inVotes = TgetInPoints(boundary_ColorFlowMap, sideCut, Tesh); % %[0.20,0.25,0.30,050,0.15]
            else
                IndxSmallFlowS1 = [IndxSmallFlowS1;frame];
                inVotes = TgetInPoints(boundary_ColorMap, sideCut, Tesh);
            end

            FInOut{frame} = inVotes > 4;   
        end
        output = FInOut;
        IndxSmallMotionFrames = IndxSmallFlowS1;
    else
        IndxSmallFlowS1 = [];
        for frame = 1:frames
            IndxImg = uint8(Images{frame});
            Indxmagflow = magflow{frame};

            Indxflow    = flow{frame};
            IndxOrgflow = Orgflow{frame};

            boundary_ColorMap     = detect_motionboundaries(IndxImg);      % not accurate
            boundary_ColorFlowMap = detect_motionboundaries(IndxImg, Indxflow);
            boundary_ColorOrgFlowMap = detect_motionboundaries(IndxImg, IndxOrgflow);

            if max(Indxmagflow(:)) > 1.0
                inVotes = TgetInPoints(boundary_ColorFlowMap, sideCut, Tesh); %[0.20,0.25,0.30,050,0.15]
                inVotesOrg = TgetInPoints(boundary_ColorOrgFlowMap, sideCut, Tesh);
            else
                IndxSmallFlowS1 = [IndxSmallFlowS1;frame];
                inVotesOrg = TgetInPoints(boundary_ColorMap, sideCut, Tesh);
                inVotes = inVotesOrg;
            end

            FInOut{frame} = inVotes > 4;   
            OrgFInOut{frame} = inVotesOrg > 4;
        end
        output = cat(3,FInOut,OrgFInOut);
        IndxSmallMotionFrames = IndxSmallFlowS1;
    end    
end
