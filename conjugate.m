% conjugate.m
% James Yates - 03/11/2019
% This script determines the conjugate of a 2D complex vector
function CONJUGATE = conjugate(z)
CONJUGATE = real(z) - 1i*imag(z);
end