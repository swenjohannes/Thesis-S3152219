function res = Lambda (H , N )
M = 2* N - 2;
C = zeros (1 , M );
G = 2* H ;
fbc = @ ( n )(( n +1).^ G + abs(n -1).^ G - 2* n .^ G )/2;
C (1: N ) = fbc (0:( N -1));
C ( N +1: M ) = fliplr ( C (2:( N -1)));
res = real ( fft ( C )).^0.5;