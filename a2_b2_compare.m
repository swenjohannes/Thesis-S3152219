%Script to create plot ..

%Parameters
parameter_set3 %Load parameter set 2
npath = 1e4; 
steps = 1200;
N = 100;
T = 0.5;

%MC simulation
[S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
[S_r, V_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);

[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 10e-4);     %obtain a2 and b2!

C0 = eta ^ 2 * (1 - exp(- kappa * T)) / ( 4 * kappa); %Compute the constant
d = 4 * kappa * theta / eta ^ 2;                      %Degrees of freedom
lambda = 4 * kappa * exp(- kappa * T) * v0 ...        %Non-centrality ..
        / (eta ^ 2 * (1 - exp(- kappa * T)));         %parameter


V_h = V_c(1201, : );                                  %Obtain final time
V_rh = V_r(1201, : );                                 %of simulated variance's

dx = 0.015;
x = (0:dx:1.5)';
ncx2 = ncx2pdf(x / C0,d,lambda);

%Create plot
figure;                                             
plot(x,ncx2,'b-','LineWidth',2); hold on;
[p,x] = hist(V_h, x); plot(x,p/sum(p));             %Empirical pdf Heston
[p,x] = hist(V_rh, x); plot(x,p/sum(p));            %Empirical pdf rough Heston
x1  = xline(a2,'--',{'a2 Heston'});
x1.LabelVerticalAlignment = 'middle';
x1.LabelHorizontalAlignment = 'center';
x2  = xline(b2,'--',{'b2 Heston'});
x2.LabelVerticalAlignment = 'middle';
x2.LabelHorizontalAlignment = 'center';
legend('Non-central chi-squared distribution','Simulated Heston values','Simulated Rough Heston values')

%Determine a2 and b2 of rough Heston
V_rh = V_r(1201, : );                                
V_rh = sort(V_rh);
a2 = V_rh(1e2);         %0.0138 param set1
b2 = V_rh(1e5 - 1e2);   %0.2676 param set1


[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants

dx = 0.001;
x = (-3:dx:3)';

figure;
l_r = log(S_c(1201, :) / S0);
l_rh = log(S_r(1201, :) / S0);
[p,x] = hist(l_r); plot(x,p/sum(p)); hold on;   %Empirical pdf rough Heston
[p,x] = hist(l_rh); plot(x,p/sum(p));    %Empirical pdf Heston

x1  = xline(a1,'--',{'a1 Heston'});
x1.LabelVerticalAlignment = 'middle';
x1.LabelHorizontalAlignment = 'center';
x2  = xline(b1,'--',{'b1 Heston'});
x2.LabelVerticalAlignment = 'middle';
x2.LabelHorizontalAlignment = 'center';
legend('Heston empirical dist','Rough Heston empirical dist')