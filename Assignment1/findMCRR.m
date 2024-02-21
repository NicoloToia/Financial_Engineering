function M = findMCRR (optionPriceBlack, F0, K, B, TTM, sigma, flag, bidAsk)
%function to find the lowest M such that error in CRR is less than Bid Ask spread
%INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% TTM:     time-to-maturity
% sigma: volatility
% flag:  1 call, -1 put
%bidAsk:  bid ask spread

M = 1; %fix M as the lowest value possible assumed by this factor in order not to lose any case
OptionPriceCRR = EuropeanOptionCRR2(F0,K,B,TTM,sigma,M,flag); 
Error = abs(OptionPriceCRR-optionPriceBlack);     
while Error>bidAsk
    M = M+1; 
    OptionPriceCRR = EuropeanOptionCRR2(F0,K,B,TTM,sigma,M,flag); 
    Error = abs(OptionPriceCRR-optionPriceBlack);     
end 
