function nDirection = gzTellDirection(aveU, aveV)
%           *
%        3  * 2
%***********************
%        4  * 1
%           *
theta    = atan2(aveV,aveU);
A_degree = theta*180/pi;

nDirection = 1;
%%
% if A_degree >= 0  & A_degree <= 45
%     nDirection =  12;
% elseif  A_degree > 45 & A_degree <= 90
%     nDirection =  23;    
% elseif  A_degree > 90 & A_degree <= 135
%     nDirection =  34;    
% elseif  A_degree > 135 & A_degree <=180
%     nDirection =  45;
% elseif  A_degree >= -180 & A_degree <= -135
%     nDirection = 56;
% elseif  A_degree > -135 & A_degree <= -90
%     nDirection = 67;
% elseif  A_degree > -90 & A_degree <= -45
%     nDirection =  78;
% elseif  A_degree > -45 & A_degree < 0
%     nDirection =  18;
% end 
%%
if A_degree >= -5  & A_degree <= 5
    nDirection =  1;
elseif  A_degree > 5 & A_degree <= 40
    nDirection = 12;    
elseif  A_degree > 40 & A_degree <= 50
    nDirection =  2;    
elseif  A_degree > 50 & A_degree <= 85
    nDirection =  23;
elseif  A_degree > 85 & A_degree <= 95
    nDirection = 3;
elseif  A_degree > 95 & A_degree <= 130
    nDirection = 34;
elseif  A_degree > 130 & A_degree <= 140
    nDirection =  4;
elseif  A_degree > 140 & A_degree <= 175
    nDirection =  45;
elseif A_degree > 175
    nDirection =  5;
elseif A_degree >= -180 & A_degree < -175
    nDirection =  5;
elseif A_degree >= -175 & A_degree < -140
    nDirection = 56;
elseif A_degree >= -140 & A_degree < -130
    nDirection =  6;
elseif A_degree >= -130 & A_degree < -95
    nDirection = 67;
elseif A_degree >= -95 & A_degree < -85
    nDirection = 7;
elseif A_degree >= -85 & A_degree < -50
    nDirection = 78;
elseif A_degree >= -50 & A_degree < -40
    nDirection =  8;
elseif A_degree >= -40 & A_degree < -5
    nDirection = 81;
end 

end 