function [a, b] = cos_truncation_range_v2( c1, c2, c4, L )

%{
 This code computes the COS Method truncation range (integration bounds) for a given model with a given L
 based as a function of its cumulants, following Fang and Oosterlee (2008). 

 We take the absolute value of c2 because, for some set of
 parameters not satisfying the Feller condition i.e.
 [2 * u_bar * lambda >= eta ^2], c2 may be negative.

 Since c4 is particularly involved for the Heston model, bounds should be based on c1 and c2 only, as
 in Fang and Oosterlee (2008).

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 [a, b] = cos_truncation_range_v2( c1, c2, c4, L )

 Inputs : c1            - first cumulant
        : c2            - second cumulant
        : c4            - 4th cumulant (if it exists)
        : L             - cumulant range scaling, inputs are 10 or 12


Outputs : a             - lower COS truncation point (lower integration bound)
        : b             - upper COS truncation point (upper integration bound)

%}

% Compute integration bounds using c1, c2 and c4 if given
a = c1 - L * sqrt( abs( c2 ) + sqrt( abs(c4) ) );
b = c1 + L * sqrt( abs( c2 ) + sqrt( abs(c4) ) );

if ( L ~= 10 && L ~= 12 )
    warning('L is recommended to be 10 for CGMY, VG and BS model. 12 is recommended for the Heston model.')
end