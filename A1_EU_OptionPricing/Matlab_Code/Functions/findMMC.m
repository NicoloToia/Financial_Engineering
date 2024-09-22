function M_MC = findMMC (optionPriceBlack, F0, K, B, TTM, sigma, flag, bidAsk)

M = 1; 
g = randn(1, 1); 
EuropeanOptionMC(F0,K,B,T,sigma,N,flag)
Error = B/sqrt(M)*std(max((Ftt - K), 0));  
while Error>bidAsk
    M = M+1; 
    OptionPriceCRR = EuropeanOptionCRR2(F0,K,B,TTM,sigma,M,flag); 
    Error = abs(OptionPriceCRR-optionPriceBlack);     
end 


g = randn(M(i),1);
    % Monte Carlo simulation (one time step)
    % Compute the value of the forward at time T for each simulation
    % Black Model: Ft = F0 * exp(-(sigma^2)*T*0.5 + sigma*sqrt(T)*g)
    Ftt = F0 * exp( -0.5 * sigma^2 * T  + sigma * sqrt(T) * g);
    
    % Compute the call option price for each simulation
    stdEstim(i)=B^2*std(max((Ftt - K),0))/sqrt(M(i));