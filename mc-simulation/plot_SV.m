function plot_SV(S, V, name, S_path, V_path)
% Uses the discritation scheme as presented in Emre Erkan 2020.
%
%  Usage:      [S, V] = rough_heston_full_mc(..);
%               
%  Inputs:      S       matrix containing simulated prices
%               V       matrix containing simulated volatilities
%               name    main title of the plot
%               S_path  range of paths to be displayed
%               V_path  range of paths to be displayed
%  Output:      a displayed plot


%Default: use all generated paths!
if (~ exist ( 'S_path' , 'var' )) S_path = 1:size(S, 2); end 
if (~ exist ( 'V_path' , 'var' )) V_path = 1:size(S, 2); end

%Plotting!
f = figure;
subplot(1,2, 1); plot(S(:, S_path)); title('Simulated asset price'); xlabel('number of steps'); ylabel('S') 
subplot(1,2, 2); plot(V(:, V_path)); title('Simulated Volatility'); xlabel('number of steps'); ylabel('V') 
sgtitle(name)

end