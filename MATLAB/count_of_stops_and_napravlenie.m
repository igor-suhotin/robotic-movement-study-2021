function [count, napr] = count_of_stops_and_napravlenie(current_k, current_myu, b)
    global k myu tol k1 k2 Vtol;
    k = current_k;
    k1=b*k;
    k2=k;
    
    myu = current_myu;
    fprintf('count_of_stops_and_napravlenie  b=%f  k=%f  myu=%f \n', b, k, myu);
    tol = 1e-14; % 0.00000001;
    Xtol=0.000000001;
    Vtol=0.00000000001; %0.0000000001; убрал 3 ноля
    x0 = 0;
    v0 = 0;
    t_start = 0;
    t_fin = 4 * pi;
    %t_step = (t_fin - t_start) / 10000; %0.01;
    t_step = 0.01;
    t = t_start:t_step:t_fin;
    out = [x0, v0];
    opts = odeset('RelTol',tol*100,'AbsTol',tol);
    
    [t, out] = ode45(@MY_equation, t, out,opts);
    
    X = out(:, 1);
    V = out(:, 2);
    
    
    if max(abs(V)) < Vtol
        count = 4; % Покой нам только снился.
    else
        [isPositiveV, isNegativeV] = deal(false); 
        count = 0;
        leng_v = length(V);
%         beginV = floor(leng_v / 2) + 1;
        i = 1;
        while i < leng_v
            isPositiveV = isPositiveV || V(i) > Vtol;
            isNegativeV = isNegativeV || V(i) < -Vtol;
            if abs(V(i)) < Vtol && abs(V(i + 1)) < Vtol 
                count = count + 1;
                while i < length(V) - 1 && abs(V(i)) < Vtol
                    i = i + 1;
                end
            end
            i = i + 1;
        end
        %if (abs(V(beginV)) < Vtol) % && (abs(V(beginV)+1) < Vtol) % && (abs(V(leng_v)) < Vtol)
%         if (abs(V(leng_v-1)) < Vtol) && (abs(V(leng_v)) < Vtol)
%             count = count - 1;
%         end
%         fprintf('count=%i  ', count);
        count = floor(count/2);
        if count == 1 && not (isNegativeV && isPositiveV)
            count = 3;
        end
    end
    leng_x = length(X);
    napr = sign(X(leng_x) - X(floor(leng_x / 2)));
    if abs(X(leng_x) - X(floor(leng_x / 2))) < Xtol
        napr = 0;
    end
%     fprintf('itog count=%i\n', count);
end