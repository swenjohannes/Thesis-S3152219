function B2 = B2_v2(w1, w2, kappa, rho, eta, tau, a1, b1, a2, b2)

w1 = w1' * pi / (b1 - a1);  %Turn w1 in w1* and make a column vector
w2 = w2  * pi / (b2 - a2);

eta2 = eta ^ 2; %Notation!
beta = kappa - 1i * rho * eta * w1; 
D = sqrt(beta .^ 2 + eta2 * w1 .* (1i + w1));

h = (beta - D - 1i * w2 * eta2) ./ (beta + D - 1i * w2 * eta2);
hemdt = h  .* exp (-D * tau); %To ease up calculations

B2 = (beta - D - (beta + D) .* hemdt) ./ (eta2 * (1 - hemdt)); %Get B2
end
 
% %Test
% w1 = 100 * pbma;
% w2 = -100 * pbma;
% tau = 0.5;