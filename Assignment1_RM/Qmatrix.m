function [Q] = Qmatrix(IG_h,HY_h)
% QMATRIX  Derive the market-implied rating transition matrix based on the
% available market data (hazard rates)
%
% INPUT
% IG_h: hazard rates for the investment grade bonds
% HY_h: hazard rates for the high yield bonds
%
% OUTPUT
% Q: market-implied rating transition matrix

% extract the hazard rates from the input
h_IG_1=IG_h(1,2);
h_IG_2=IG_h(2,2);
h_HY_1=HY_h(1,2);
h_HY_2=HY_h(2,2);
% compute the one year transition probabilities
q13 = 1 - exp(-h_IG_1);
q23 = 1 - exp(-h_HY_1);
% compute the two year transition probabilities
q13_2 = 1 - exp(-h_IG_2-h_IG_1);
q23_2 = 1 - exp(-h_HY_2-h_HY_1);
% solve the system of equations to find the market-implied transition matrix
A =[ q13 q23 0 0; 
    0 0 q13 q23;
    1 1 0 0;
    0 0 1 1 ];
b = [q13_2-q13 q23_2-q23 1-q13 1-q23]';
x = A\b;
% construct the market-implied transition matrix
Q = [ x(1) x(2) q13;
    x(3) x(4) q23;
    0 0 1];

end 

