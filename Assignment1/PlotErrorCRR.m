function [nCRR,errCRR]=PlotErrorCRR(F0,K,B,T,sigma)
% error plot for CRR method
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility


% calcolo iterazioni ottimali CRR ed errore
clc
clear all
close all
errCRR=1;
nCRR=100;
Eta=10^-3;

while errCRR>Eta
    for i=1:2
        pricingMode = i; % 1 ClosedFormula, 2 CRR
        OptionPrice(i) = EuropeanOptionPrice(F0,K,B,T,sigma,pricingMode,nCRR);
    end
    errCRR=abs(OptionPrice(1)-OptionPrice(2));
    nCRR=nCRR+5;
end

N=[0:5:nCRR];
for j=1:length(N)
    for i=1:2
        pricingMode = i; % 1 ClosedFormula, 2 CRR
        OptionPrice(i) = EuropeanOptionPrice(F0,K,B,T,sigma,pricingMode,nCRR);
    end
    error(j)=abs(OptionPrice(1)-OptionPrice(2));
end
loglog(N,error)

