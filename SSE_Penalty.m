function [SSE] = SSE_Penalty(c, y, lambda, basis, dt)

% Calculate the least squared error with a roughness penalty for chi and xi
% Inputs: 
%   c: a column vector of coefficients
%   y: a column vector of function values (chi or xi)
%   lambda: smoothing parameter
%   basis: the basis for functions chi(t) and xi(t)
%   dt: delta t
% Outputs:
%   SSE: sum squared error with penalty

[N, n] = size(knots);

D2basis = zeros(N-2, n); % second order derivatives

for i = 2: N-1
    D2basis(i-1, :) = (basis(i+1, :) - 2*basis(i, :) + basis(i-1, :))./ dt^2 ;
end

R = D2basis' * D2basis * dt; 

J = c' * R * c; 
SSE = (y - basis * c)' * (y - basis * c) + lambda * J;



