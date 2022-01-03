function [c,ceq] = confuneq(x) 
% Nonlinear inequality constraints
c = 2 * x(2) * x(3) - x(4)^2; %FELLER CONDITION
% Nonlinear equality constraints
ceq = [];
end