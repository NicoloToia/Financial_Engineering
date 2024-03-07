function [Q] = matrix(R, IG_z_curve, HY_z_curve)

%    Derive the market-implied rating transition matrix based on the
%   available market data (ZC risk-free curve and risky bond prices)

q13 = (1 - exp(-IG_z_curve(1))) / (1 - R);
q23 = (1 - exp(-HY_z_curve(1))) / (1 - R);

q13_2 = ( 1 - exp(-IG_z_curve(2) - IG_z_curve(1)) ) / (1 - R);
q23_2 = ( 1 - exp(-HY_z_curve(2) - HY_z_curve(1)) ) / (1 - R);

A=[q13 q23 0 0; 
   0 0 q13 q23;
   1 1 0 0;
   0 0 1 1];

b= [q13_2-q13 q23_2-q23 1-q13 1-q23]';

x=A\b;

Q=[x(1) x(2) q13; x(3) x(4) q23; 0 0 1];
end 

