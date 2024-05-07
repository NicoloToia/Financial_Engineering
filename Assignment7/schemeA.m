function [p_u, p_m, p_d] = schemeA(l, mu)
% Scheme A: central scheme of the Trinomial Tree for interest rates
%
% INPUTS:
%   l: index of the node
%   mu: drift of the interest rate process
%
% OUTPUTS:
%   p_u: probability of the interest rate going up
%   p_m: probability of the interest rate staying the same
%   p_d: probability of the interest rate going down

% compute the probabilities of the trinomial tree
p_u = 1/2*(1/3 - mu*l + mu^2*l.^2);
p_m = 2/3 - mu^2*l^2;
p_d = 1/2*(1/3 + mu*l + mu^2*l^2);

end