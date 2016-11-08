function iou = inters_union(box1,box2)

% AREA = RECTINT(A,B) returns the area of intersection of the rectangles specified by position vectors A and B.  
% 'A' position vector: is a four-element vector [X,Y,WIDTH,HEIGHT],
% where the point defined by X and Y specifies one corner of the rectangle, and WIDTH and HEIGHT define the size in units along the x-and y-axes respectively.

inters = rectint(box1,box2); %returns the area of intersection of the rectangles specified by position vectors box1 and box2.
ar1 = box1(:,3).*box1(:,4);
ar2 = box2(:,3).*box2(:,4);
union = bsxfun(@plus,ar1,ar2')-inters;

iou = inters./(union+eps);
