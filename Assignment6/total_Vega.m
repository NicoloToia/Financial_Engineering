function Vega = total_Vega(SpotVols)

    % shock the volatility of the caps by 1bp
    bp = 0.0001;
    shock = bp * ones(size(SpotVols));
    shockVols = SpotVols + shock;

    % compute the total vega
    Vega = 

end