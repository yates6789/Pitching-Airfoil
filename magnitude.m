% magnitude.m
% James Yates - 03/11/2019
% This script determines the magnitude of a 2D complex vector
function MAGNITUDE = magnitude(z)
MAGNITUDE = sqrt((real(z).*real(z)) + (imag(z).*imag(z)));
end