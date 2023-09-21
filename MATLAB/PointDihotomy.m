% function m = PointDihotomy( b, k, myu_begin, C_begin, myu_end, C_end, steps)
function P = PointDihotomy( b, k_begin, myu_begin, C_begin, k_end, myu_end, C_end, steps)

%     global k myu tol k1 k2 Vtol;
    
    k=(k_begin+k_end)/2;
    myu = (myu_begin+myu_end)/2;
    
    if steps > 0
        C = count_of_stops_and_napravlenie(k, myu, b);
        if (C(1)==C_begin(1))
            P=PointDihotomy(b, k, myu, C, k_end, myu_end, C_end, steps-1);
%             m=PointDihotomy( b, k, myu, C, myu_end, C_end, steps-1);
        else
%             m=PointDihotomy( b, k, myu_begin, C_begin, myu, C, steps-1);
            P=PointDihotomy(b, k_begin, myu_begin, C_begin, k, myu, C, steps-1);
        end
    else
        P=[k,myu];
        fprintf('Dihotomy  k=%f  myu=%f', k, myu);
    end
        
end

