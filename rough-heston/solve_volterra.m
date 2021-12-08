function [y,Dalpha_y] = SolveVIE(f, y0, alpha, T, N)

%{
    Description: Solves the Volterra integral equation (VIE)         
                  y(t) = (1/gamma(alpha))*int_0^t(t-u)^(alpha - 1)f(u,y(u))du
    
                  The equation is solved on the equidistant grid
                  0 < T/N < 2*T/N < 3*T/N < ... < N*T/N = T
    
                  The notation is the same as used in (Diethelm, 2004).
    
    Parameters:
      f:      [function] Function of two variables from the VIE. We allow it 
              to return a [Mx1] vector but then it should also accept [Mx1]
              vectors for its second argument.
      y0:     [1x1 real] Initial value 
      alpha:  [1x1 real] Number between 0 and 1.
      T:      [1x1 real] Upper time point to solve VIE on.
      N:      [1x1 integer] Number steps to discretize [0,T] into.
    
    Output: 
      y:        [Mx(N+1) real] Solution of VIE, 
                i.e. [y(0),y(T/N),y(2*T/N),...,y(T)]
      Dalpha_y: [Mx(N+1) real] Fractional (alpha) derivative of solution, 
                i.e. [D^(alpha)y(0),D^(alpha)y(T/N),...,D^(alpha)y(T)]
    
    References:
      - Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error
      analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52,
      2004.
%}

% Initialization:
h = T / N;
t = (0:h:T);
[y,Dalpha_y] = deal(NaN(M,N+1));


% Run scheme:
y(:,1) = y0;
Dalpha_y(:,1) = f(0,y0);
for k=0:N-1
    js = (0:1:k);

    % Compute predictor:
    yp = sum(b_j_kp1(js,k).*Dalpha_y(:,1:k+1),2)./gamma(alpha);

    % Compute solution:
    if k==0
        y(:,2) = (a_0_kp1(k)*Dalpha_y(:,1) ...
            + a_kp1_kp1*f(t(k+2),yp))./gamma(alpha);
    else
        y(:,k+2) = (a_0_kp1(k)*Dalpha_y(:,1) ...
            + sum(Dalpha_y(:,2:k+1).*a_j_kp1(js(2:end),k),2)...
            + a_kp1_kp1*f(t(k+2),yp))./gamma(alpha);
    end

    % Compute fractional derivative:
    Dalpha_y(:,k+2) = f(h,y(:,k+2));

end

if any(any(isnan(y))) || any(any(isnan(Dalpha_y)))
    error('SolveVIE: NaNs produced!');
end

end


%% Internal functions 
function a = a_j_kp1(j, k)
%{
    Description:  Computes the constant a_j_kp1 

    Parameters:
      j:      
      k:     [1x1 real] 
      alpha:  [1x1 real] Number between 0 and 1.

    Output: 
      a:        [1x(k+1)] vector of a_j_kp1's
                i.e. [y(0),y(T/N),y(2*T/N),...,y(T)]
%}

% Define coefficient functions:
term1 = ( (h^alpha) / (alpha*(alpha + 1)) );


a_0_kp1 = @(k)( dummy1 *( k.^(alpha+1)-(k - alpha).*((k+1).^(alpha))));
a_j_kp1 = @(j,k)( dummy1*( (k-j+2).^(alpha+1) ...
    + (k-j).^(alpha+1) - 2*(k-j+1).^(alpha+1) ) );
a_kp1_kp1 = dummy1;
end

function b = b_j_kp1(j, k, h, alpha)

%{
    Description:  Computes the constant a_j_kp1 

    Parameters:
      j:      [1xN real] vector of j's to compute the function for
      k:      [1x1 real] The k to compute the function for 
      h:      [1x1 real] The grid width
      alpha:  [1x1 real] Number between 0 and 1.

    Output: 
      a:        [1x(k+1)] vector of a_j_kp1's
                i.e. [y(0),y(T/N),y(2*T/N),...,y(T)]
%}
b =  (h^alpha)/alpha *  ((k+1-j).^(alpha) - (k-j).^(alpha)); 
end

