function int = integral_jam(t0,talfa,ti,sigma,a)

for i = 1:length(ti)
    ti(i) = yearfrac(t0,ti(i),3);
end

talfa = yearfrac(t0,talfa,3);

integrand = @(u) ((sigma/a * (1 - exp(-a * (ti-u)))).^2 - (sigma/a * (1 - exp(-a * (talfa-u)))).^2 );

int = integral(integrand, 0, talfa,"ArrayValued",true);

end