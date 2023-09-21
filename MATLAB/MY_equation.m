function Xdot = MY_equation(t,X)
global myu k1 k2 Vtol
n=cos(t)+myu;

s=-sin(t);

if abs(X(2)) < Vtol
    X(2) = 0;
    if s < -k1*n
        fc = -k1*n;
    elseif s > k2*n
        fc = k2*n;
    else 
        fc = s;
    end      
elseif X(2) > 0
        fc = -k1*n;
else
        fc = k2*n;
end

Xdot(1)=X(2);
Xdot(2)=fc-s;
Xdot=Xdot';
end
