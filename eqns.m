function f = eqns(p,kp,p0)
f(1) = kp(1)*(p(2)^2)- p(1);
f(2) = kp(2)*(p(4)^2)- p(3);
f(3) = (21/79)*(2*p(1)+p(2)+p(5)+p(8))-(2*p(3)+p(4)+p(5)+p(7));
f(4) = p0 - p(1) - p(2) - p(3) - p(4) - p(5) - p(6) - p(7) - p(8);
f(5) = kp(3)*(p(2)*p(4)) - p(5);
f(6) = kp(4)*(p(7)*p(6)) - p(4);
f(7) = kp(5)*(p(8)*p(6)) - p(2);
f(8) = p(6) - p(7) - p(8);


