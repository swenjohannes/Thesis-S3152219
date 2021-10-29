function res = FGN ( lambda , NT )
if (~ exist ( 'NT' , 'var' ))
    NT = 1;
end

M = size(lambda ,2);
a = bsxfun( @times , ifft (randn( NT , M ) ,[] ,2) , lambda );
res = real ( fft (a ,[] ,2));
res = res(: ,1:( M /2));
