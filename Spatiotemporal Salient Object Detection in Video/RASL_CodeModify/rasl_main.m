function [D, Do, A, E, xi, numIterOuter, numIterInner ] = rasl_main(fileNames, transformations, numImages, raslpara, destDir)

% ---------------------------------------------------
% Batch image alignment: RASL main loop
%
% input: fileNames               --- the names of images
%        transformations         --- the initial transformation matrix
%        numImages               --- the number of images
%        raslpara                --- parameters for RASL
%        destDir                 --- output
%
% output: D                      --- input images in canonical frame
%         Do                     --- output aligned images
%         A                      --- low-rank component
%         E                      --- sparse error component
%         xi                     --- transformation parameters
%         numIterOuter           --- number of outer loop iterations
%         numIterInner           --- total number of inner loop iterations
% ---------------------------------------------------

%% read and store full images

APGorALM_flag = 0;
% APGorALM_flag: 1: APG algorithm
%                0: inexact ALM algorithm  

fixGammaType = 1 ;

if ~fixGammaType
    if exist(fullfile(rootPath, userName, 'gamma_is_ntsc'), 'file')
        gammaType = 'ntsc' ;
    elseif exist(fullfile(rootPath, userName, 'gamma_is_srgb'), 'file')
        gammaType = 'srgb' ;
    elseif exist(fullfile(rootPath, userName, 'gamma_is_linear'), 'file')
        gammaType = 'linear' ;
    else
        error('Gamma type not specified for training database!  Please create a file of the form gamma_is_*') ;
    end
else
    gammaType = 'linear' ;
end

sigma0 = 2/5 ;
sigmaS = 1 ;

deGammaTraining = true ;

I0 = cell(raslpara.numScales,numImages) ; % images
I0x = cell(raslpara.numScales,numImages) ; % image derivatives
I0y = cell(raslpara.numScales,numImages) ;

for fileIndex = 1 : numImages
    
    currentImage = double(imread(fileNames{fileIndex}));
    
    % Use only the green channel in case of color images    
    if size(currentImage,3) > 1,   currentImage = currentImage(:,:,2);            end
    if deGammaTraining,      currentImage = gamma_decompress(currentImage, gammaType); end
    
    currentImagePyramid = gauss_pyramid( currentImage, raslpara.numScales,...
        sqrt(det(transformations{fileIndex}(1:2,1:2)))*sigma0, sigmaS );
        
    for scaleIndex = raslpara.numScales:-1:1
        I0{scaleIndex,fileIndex} = currentImagePyramid{scaleIndex};
        
        % image derivatives
        I0_smooth = I0{scaleIndex,fileIndex};
        I0x{scaleIndex,fileIndex} = imfilter( I0_smooth, (-fspecial('sobel')') / 8 );
        I0y{scaleIndex,fileIndex} = imfilter( I0_smooth,  -fspecial('sobel')   / 8 );
    end   
end



%% get the initial input images in canonical frame

imgSize = raslpara.canonicalImageSize ; 

xi_initial = cell(1,numImages) ; % initial transformation parameters
for i = 1 : numImages
    if size(transformations{i},1) < 3
        transformations{i} = [transformations{i} ; 0 0 1] ;
    end
    xi_initial{i} = projective_matrix_to_parameters(raslpara.transformType,transformations{i});
end

D = [] ;
for fileIndex = 1 : numImages
    % transformed image
    Tfm = fliptform(maketform('projective',transformations{fileIndex}'));
            
    I   = vec(imtransform(I0{1,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    y   = I; 
    y = y / norm(y) ; % normalize

    D = [D y] ;
end

if raslpara.saveStart
    save(fullfile(destDir, 'original.mat'),'D','xi_initial') ;
end


%% start the main loop

frOrig = cell(1,numImages) ;
T_in = cell(1,numImages) ;

T_ds = [ 0.5,   0, -0.5; ...
         0,   0.5, -0.5   ];
T_ds_hom = [ T_ds; [ 0 0 1 ]];

numIterOuter = 0 ; 
numIterInner = 0 ;

tic % time counting start

for scaleIndex = raslpara.numScales:-1:1 % multiscale
    
    iterNum = 0 ;  % iteration number of outer loop in each scale
    converged = 0 ;
    prevObj = inf ; % previous objective function value
    
    imgSize = raslpara.canonicalImageSize / 2^(scaleIndex-1) ;    
    xi = cell(1,numImages) ;
    
    for fileIndex = 1 : numImages
             
        if scaleIndex == raslpara.numScales
            T_in{fileIndex} = T_ds_hom^(scaleIndex-1)*transformations{fileIndex}*inv(T_ds_hom^(scaleIndex-1)) ;
        else
            T_in{fileIndex} = inv(T_ds_hom)*T_in{fileIndex}*T_ds_hom ;
        end
        
        % for display purposes
        if raslpara.DISPLAY > 0
            fr = [1 1          imgSize(2) imgSize(2) 1; ...
                  1 imgSize(1) imgSize(1) 1          1; ...
                  1 1          1          1          1 ];
            
            frOrig{fileIndex} = T_in{fileIndex} * fr;
        end
        
    end
    
    while ~converged

        iterNum = iterNum + 1 ;
        numIterOuter = numIterOuter + 1 ;
        
        D = [] ; J = cell(1,numImages) ;
        disp(['Scale ' num2str(scaleIndex) '  Iter ' num2str(iterNum)]) ;
        
        for fileIndex = 1 : numImages

            % transformed image and derivatives with respect to affine parameters
            Tfm = fliptform(maketform('projective',T_in{fileIndex}'));
            
            I   = vec(imtransform(I0{scaleIndex,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            Iu  = vec(imtransform(I0x{scaleIndex,fileIndex},Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            Iv  = vec(imtransform(I0y{scaleIndex,fileIndex},Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            y   = I; %vec(I);

            Iu = (1/norm(y))*Iu - ( (y'*Iu)/(norm(y))^3 )*y ;
            Iv = (1/norm(y))*Iv - ( (y'*Iv)/(norm(y))^3 )*y ;

            y = y / norm(y) ; % normalize
            D = [D y] ;

            % transformation matrix to parameters
            xi{fileIndex} = projective_matrix_to_parameters(raslpara.transformType,T_in{fileIndex}) ; 
            
            % Compute Jacobian
            J{fileIndex} = image_Jaco(Iu, Iv, imgSize, raslpara.transformType, xi{fileIndex});
        end
        

        lambda = raslpara.lambdac/sqrt(size(D,1)) ; 

        
        % RASL inner loop
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % using QR to orthogonalize the Jacobian matrix
        for fileIndex = 1 : numImages
            [Q{fileIndex}, R{fileIndex}] = qr(J{fileIndex},0) ;
        end
        
        if APGorALM_flag == 1
            [A, E, delta_xi, numIterInnerEach] = rasl_inner_apg(D, Q, lambda, raslpara.inner_tol, raslpara.inner_maxIter, raslpara.continuationFlag, raslpara.mu) ;
        else
            [A, E, delta_xi, numIterInnerEach] = rasl_inner_ialm(D, Q, lambda, raslpara.inner_tol, raslpara.inner_maxIter);
        end
        for fileIndex = 1 : numImages
            delta_xi{fileIndex} = inv(R{fileIndex})*delta_xi{fileIndex} ;
        end
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        

        numIterInner = numIterInner + numIterInnerEach ;

        curObj = norm(svd(A),1) + lambda*norm(E(:),1) ;
        disp(['previous objective function: ' num2str(prevObj) ]);
        disp([' current objective function: ' num2str(curObj) ]);

        % step in paramters
        for i = 1 : numImages
            xi{i} = xi{i} + delta_xi{i};
            T_in{i} = parameters_to_projective_matrix(raslpara.transformType,xi{i});
        end
        
        % save intermedia results
        if raslpara.saveIntermedia
            matName = strcat('scale_', num2str(scaleIndex),'_iter_', num2str(iterNum),'.mat') ;
            save(fullfile(destDir, matName),'D','A','E','xi') ;
        end

        if raslpara.DISPLAY > 0
            for i = 1 : numImages
                figure(1); clf ;
                imshow(I0{scaleIndex,i},[],'Border','tight');
                hold on;
                
                Tfm = fliptform(maketform('projective',inv(T_in{i}')));
                curFrame = tformfwd(fr(1:2,:)', Tfm )';
                plot( frOrig{i}(1,:),   frOrig{i}(2,:),   'g-', 'LineWidth', 2 );
                plot( curFrame(1,:), curFrame(2,:), 'r-', 'LineWidth', 2 );
%                 hold off;
%                 print('-f1', '-dbmp', fullfile(destDir, num2str(i))) ;
               
            end
        end
        
        if ( (prevObj - curObj < raslpara.stoppingDelta) || iterNum >= raslpara.maxIter )
            converged = 1;
            if ( prevObj - curObj >= raslpara.stoppingDelta )
                disp('Maximum iterations reached') ;
            end
        else
            prevObj = curObj;
        end
        
    end
end

timeConsumed = toc

disp(['total number of iterations: ' num2str(numIterInner) ]);
disp(['number of outer loop: ' num2str(numIterOuter) ]);

%% save the alignment results

Do = [] ;
for fileIndex = 1 : numImages

    Tfm = fliptform(maketform('projective',T_in{fileIndex}'));
            
    I   = vec(imtransform(I0{1,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    y   = I; 
    y = y / norm(y) ; % normalize

    Do = [Do y] ;
end

if raslpara.saveEnd
    save(fullfile(destDir, 'final.mat'),'Do','A','E','xi') ;
end


outputFileName = fullfile(destDir, 'results.txt'); 
fid = fopen(outputFileName,'a') ;
fprintf(fid, '%s\n', [' total number of iterations: ' num2str(numIterInner) ]) ;
fprintf(fid, '%s\n', [' number of outer loop ' num2str(numIterOuter) ]) ;
fprintf(fid, '%s\n', [' consumed time: ' num2str(timeConsumed)]) ;
fprintf(fid, '%s\n', [' the parameters :']) ;
fprintf(fid, '%s\n', [' transformType ' raslpara.transformType ]) ;
fprintf(fid, '%s\n', [' lambda ' num2str(raslpara.lambdac) ' times sqrt(m)']) ;
fprintf(fid, '%s\n', [' stoppingDelta of outer loop ' num2str(raslpara.stoppingDelta) ]) ;
fprintf(fid, '%s\n', [' stoppingDelta of inner loop ' num2str(raslpara.inner_tol)]) ;
if APGorALM_flag == 1
    fprintf(fid, '%s\n', [' optimization in inner loop is using APG algorithm']) ;
    fprintf(fid, '%s\n', [' continuationFlag  ' num2str(raslpara.continuationFlag) ]) ;
    fprintf(fid, '%s\n', [' mu of inner loop ' num2str(raslpara.mu) ]) ;
else
    fprintf(fid, '%s\n', [' optimization in inner loop is using inexact ALM algorithm']) ;
end
fclose(fid);

