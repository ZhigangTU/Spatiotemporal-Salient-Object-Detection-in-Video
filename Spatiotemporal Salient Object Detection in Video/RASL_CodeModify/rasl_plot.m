function rasl_plot(destDir, numImage, canonicalImageSize, layout)

%% load in data
% initial input images
load(fullfile(destDir, 'original.mat'), 'D') ;

% alignment results
load(fullfile(destDir, 'final.mat'), 'Do','A','E') ;

%% display

% layout
if nargin < 4
    xI = ceil(sqrt(numImage)) ;
    yI = ceil(numImage/xI) ;

    gap = 2;
    gap2 = 1; % gap2 = gap/2;
else
    xI = layout.xI ;
    yI = layout.yI ;

    gap = layout.gap ;
    gap2 = layout.gap2 ; % gap2 = gap/2;
end
container = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap); 
% white edges
bigpic = cell(xI,yI); % (xI*canonicalImageSize(1),yI*canonicalImageSize(2));

% D
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(D(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(D))],'Border','tight')
title('Input images') ;

% Do
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(Do(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(Do))],'Border','tight')
title('Aligned images') ;


% A
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(A(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(A))],'Border','tight')
title('Aligned images adjusted for sparse errors') ;

% E
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(E(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure
imshow(abs(cell2mat(bigpic)),[],'DisplayRange',[0 max(max(abs(E)))],'Border','tight')
title('Sparse corruptions in the aligned images') ;

% average face of D, Do and A
% bigpic = cell(1,3); 
% container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(sum(D,2), canonicalImageSize);
% bigpic{1,1} = container;
% container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(sum(Do,2), canonicalImageSize);
% bigpic{1,2} = container;
% container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(sum(A,2), canonicalImageSize);
% bigpic{1,3} = container;
% 
figure
subplot(1,3,1)
imshow(reshape(sum(D,2), canonicalImageSize),[])
title('average of unaligned D')
subplot(1,3,2)
imshow(reshape(sum(Do,2), canonicalImageSize),[])
title('average of aligned D')
subplot(1,3,3)
imshow(reshape(sum(A,2), canonicalImageSize),[])
title('average of A')