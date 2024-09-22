function Caplet = CapletBachelier(fwd_Libor, Strike, Vol, T_payment, T_exercise, t0, DF_payment)
% CapletBachelier: Compute the price of a caplet using the Bachelier formula
%
% INPUT
%   fwd_Libor : Forward Libor rate between T_exercise and T_payment
%   Strike : Strike price of the caplet as a percentage change (eg. 5 for 5%)
%   Vol : Volatility of the caplet
%   T_payment : Date in which the Euribor is paid (payment date)
%   T_exercise : Exercise date of the caplet. When:
%       - the option is exercised or not
%       - Libor will be fixed
%   t0 : Settlement date
%   DF_payment : Discount factor from t0 to T_payment

% define the year fraction conventions
ACT_360 = 2;
ACT_365 = 3;
% compute the yearfrac of the option (from t0 to T_exercise)
delta_ttm = yearfrac(t0, T_exercise, ACT_365);
% compute the yearfrac between which the Libor will accrue (from T_exercise to T_payment)
delta_libor = yearfrac(T_exercise, T_payment, ACT_360);

% compute the real strike of the caplet
K = Strike/100;

% compute the argument of the normal cdf
d_n = (fwd_Libor - K) / (Vol * sqrt(delta_ttm));

% compute the two terms of the caplet
term_1 = (fwd_Libor - K) * normcdf(d_n);
term_2 = Vol * sqrt(delta_ttm) * normpdf(d_n);

% compute the price of the caplet
Caplet = DF_payment * delta_libor * (term_1 + term_2);

end