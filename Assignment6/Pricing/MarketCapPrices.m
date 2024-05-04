function CapPrices = MarketCapPrices(ttms, strikes, Vols, caplet_ttms, caplet_yf, caplet_DF, Libor)
% MARKETCAPPRICES computes the cap prices from the market data
%
% INPUTS
%   ttms: time to maturities
%   strikes: strikes
%   Vols: market volatilities
%   caplet_ttms: caplet maturities
%   caplet_yf: caplet year fractions
%   caplet_DF: caplet discount factors
%  Libor: forward Libor rates

% save the caps
CapPrices = zeros(size(Vols));

for i = 1:length(ttms)

    % find the corresponding caplet data (remember to skip the first)
    relevant_ttms = caplet_ttms(1:4*ttms(i)-1);
    relevant_yf = caplet_yf(1:4*ttms(i)-1);
    relevant_DF = caplet_DF(1:4*ttms(i)-1);
    relevant_Libor = Libor(1:4*ttms(i)-1);

    % for each strike
    for j = 1:length(strikes)
        % apply the bachelier formula in the vector form
        Strike = strikes(j) / 100 * ones(size(relevant_ttms));
        dn = (relevant_Libor - Strike) ./ (Vols(i, j) * sqrt(relevant_ttms));

        term_1 = (relevant_Libor - Strike) .* normcdf(dn);
        term_2 = Vols(i, j) * sqrt(relevant_ttms) .* normpdf(dn);

        CapPrices(i, j) = sum(relevant_DF .* relevant_yf .* (term_1 + term_2));

    end
    
end

end