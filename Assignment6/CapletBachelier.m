function Caplet = CapletBachelier(Libor, Strike, Vol, delta, T, t0, discount)

% compute the yearfrac of the ttm
ACT_365 = 3;
delta_ttm = yearfrac(t0, T, ACT_365);

d_n = (Libor - Strike) / (Vol * sqrt(delta_ttm));

Caplet = delta * discount * ( ...
    (Libor - Strike) * normcdf(d_n) + ...
    Vol * sqrt(delta_ttm) * normpdf(d_n) ...
)

end