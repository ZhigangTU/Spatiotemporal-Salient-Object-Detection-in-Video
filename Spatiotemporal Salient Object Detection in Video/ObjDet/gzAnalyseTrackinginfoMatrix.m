function [bConsistent, MaxDiffangle] = gzAnalyseTrackinginfoMatrix(TrackinfoMatrix, Nframes)
%gzAnalyseTrackinginfoMatrix(TrackinfoMatrix);

bConsistent=0;
MaxDiffangle = 0;

[m,n] = size(TrackinfoMatrix);

TrackF = 5;    % 5, 10, 15(raincar case), 30(treeman case)
TrackN = min(TrackF, Nframes);

if m >= TrackN 
    hstart = TrackinfoMatrix(1,2);
    hend   = TrackinfoMatrix(m,2);
    
    wstart = TrackinfoMatrix(1,3);
    wend   = TrackinfoMatrix(m,3);
    
    vect_u = TrackinfoMatrix(:,4);
    vect_v = TrackinfoMatrix(:,5);
    
    vect_norm = sqrt(vect_u.^2 + vect_v.^2);
    Normsort  = sort(vect_norm,'descend');
    dThresh   = Normsort(floor(m/2));%-0.001; 
    norm_good = vect_norm>dThresh;
    
    MaxDiffangle = 0;
    for j = 1:m,
        if norm_good( j )
            vector1=[vect_u(j) vect_v(j)];
            for j1 = 1:m,
                if norm_good( j1 )
                    vector2=[vect_u(j1 ) vect_v(j1)];
                    Cosvalue=dot(vector1,vector2)/(norm(vector1)*norm(vector2));
                    theta_C=acos(Cosvalue);
                    A_degree = theta_C*180/pi;
                    if A_degree>MaxDiffangle
                        MaxDiffangle = A_degree;
                    end
                end
            end
        end
    end
     
    posU = vect_u>0;
    posUratio = sum(posU)/m;
    
    posV = vect_v>0;
    posVratio = sum(posV)/m;
%     % treeman case
%     if (abs(hend - hstart) > 4 || abs(wend - wstart)> 4) && (posUratio>0.8 || posUratio<0.2 || posVratio>0.8 || posVratio<0.2)
    % raincar case
    if (abs(hend - hstart) > 3 || abs(wend - wstart)> 3) && (posUratio>0.7 || posUratio<0.3 || posVratio>0.7 || posVratio<0.3)
        bConsistent = 1;
    end
end

end