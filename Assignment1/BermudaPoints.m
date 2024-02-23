%% Bermudan Option (Point g)

%   Price also -with the Tree- a Bermudan option, where the holder has also the right to exercise 
%   the option at the end of every month, obtaining the stock at the strike price. 
M = 1000;
% S0 = K;  % TO BE CHECKED ---->   NON CREDO SERVANO
% F0 = S0*exp(-d*TTM)/B;

OptionPriceBermudan = BermudanOptionCRR(F0, K, B, TTM, sigma, M, flag);
fprintf('\nCRRPriceBermudan      :   %.4f \n',OptionPriceBermudan);

%%  Vary the Dividend Yield (Point h)

%   Pricing the Bermudan option, vary the dividend yield between 0% and 6% and compare 
%   with the corresponding European price. Discuss the results.

d=0:0.005:0.06;
Delta_price=zeros(length(d),1);
OptionPriceBermudan = zeros(length(d),1);
OptionCRR = zeros(length(d),1);

for i=1:length(d)
    F0=S0*exp(-d(i)*TTM)/B;
    OptionPriceBermudan(i) = BermudanOptionCRR(F0, K, B, TTM, sigma, M, flag);
    OptionCRR(i) = EuropeanOptionCRR(F0, K, B, TTM, sigma, M, flag);
    Delta_price=OptionPriceBermudan-EuropeanOptionCRR(F0, K, B, TTM, sigma, M, flag);
end

if OptionPriceBermudan - OptionCRR >= 0
    disp('true')
end

figure
plot(d, OptionPriceBermudan, '-xb', 'LineWidth', 1); 
hold on; 
plot(d, OptionCRR, '-or', 'LineWidth', 1);
xlabel('Dividend Yield'); 
ylabel('Option Price'); 
title('Comparison of Bermudan Option Price and European Option Price'); 
legend('Bermudan Option Price', 'European Option Price'); 