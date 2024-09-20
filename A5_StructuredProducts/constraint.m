function [c, ceq] = constraint(x, alpha)
    % function that represents the constraint in the Least Squares Calibration
    %p1=sigma
    %p2 k, p3 eta
    c = - (1-alpha) / (x(2)*x(1)^2) - x(3);
    ceq = [];
end