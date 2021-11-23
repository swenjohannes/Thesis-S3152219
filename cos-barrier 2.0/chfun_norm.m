function y = chfun_norm(sigma, r, T, w)

% Characteristic function of BSM.
% y = chfun_norm(s0, v, r, t, w)
% Inputs:
% s0: stock price
% sigma: volatility 
% r: risk-free rate
% t: time to maturity
% w: points at which to evaluate the function

mean = (r-sigma ^ 2 / 2)*T; % mean
var = sigma ^ 2 * T; % variance
y = exp( 1i * w *  mean -(w .* w * var * .5)); % characteristic function of log (St) evaluated at points w
end
