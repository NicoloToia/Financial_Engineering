function [Q] = matrix(IG_h,HY_h)

%    Derive the market-implied rating transition matrix based on the
%   available market data (ZC risk-free curve and risky bond prices)

h_IG_1=IG_h(1,2);
h_IG_2=IG_h(2,2);

h_HY_1=HY_h(1,2);
h_HY_2=HY_h(2,2);

q13 = 1 - exp(-h_IG_1);
q23 = 1 - exp(-h_HY_1);

q13_2 = 1 - exp(-2*h_IG_2);
q23_2 = 1 - exp(-2*h_HY_2);

A=[q13 q23 0 0; 
   0 0 q13 q23;
   1 1 0 0;
   0 0 1 1];

b= [q13_2-q13 q23_2-q23 1-q13 1-q23]';

x=A\b;

Q=[x(1) x(2) q13; x(3) x(4) q23; 0 0 1];
end 

