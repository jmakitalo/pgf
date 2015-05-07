function writeComplexArray(file, array)
%WRITECOMPLEXARAY   Writes a specific type of array to a file in the
%                   Fortran notion of complex numbers.
%   file: handle for a file opened for writing.
%   array: 4-dimensional array of complex numbers.

%fprintf(file, '(%.15g, %.15g)\n', real(array), ...
%                imag(array));

for n=1:size(array,4)
    for m=1:size(array,3)
         for k=1:size(array,2)
           for l=1:size(array,1)
             fprintf(file, '(%.15g, %.15g)\n', real(array(l, k, m, n)), ...
                 imag(array(l, k, m, n)));
           end
         end
    end
end
end