function [SSE, y_hat] = SSE_Penalty_FR(c, y, chi, xi, lambda, basis1, coe1, basis2, dt)

% Calculate the least squared error with a roughness penalty for functional regression
% Inputs: 
%   c: a column vector of coefficients, [beta0, beta1, beta2, c^{alpha1}, c^{alpha2}]
%   y: a column vector of function values
%   chi: time series chi
%   xi: time series xi
%   lambda: smoothing parameter
%   basis1: B_{i,m1}^{(1)}, the basis for functions chi(t) and xi(t)
%   coe1: [c^{chi}, c^{xi}]
%   basis2: B_{j,m2}^{(2)}, the basis for functions alpha1(t) and alpha2(t)
%   dt: delta t
% Outputs:
%   SSE: sum squared error with penalty
%   y_hat: estimates of y

beta0 = c(1);
beta1 = c(2);
beta2 = c(3);

n = (length(c) - 3) / 2; % the length of c^{alpha1} and c^{alpha2}
coe_func_chi = c(4: 3+n)'; % c^{alpha1}
coe_func_xi = c(4+n: end)'; % c^{alpha2}
c_hat_chi = coe1(:, 1); % c^{chi}
c_hat_xi = coe1(:, 2); % c^{xi}

SSE = 0;
N = length(y);
y_hat = zeros(N, 1); 

func_chi = (basis1 * c_hat_chi) .* (basis2 * coe_func_chi); % chi(t)
func_xi = (basis1 * c_hat_xi) .* (basis2 * coe_func_xi); % xi(t)

for i = 1: N
    y_hat(i) = beta0 + beta1*chi(i) + beta2*xi(i) + sum(func_chi(i:end))*dt + sum(func_xi(i:end))*dt;
    SSE = SSE + (y(i) - y_hat(i)).^2;
end

D2yhat = zeros(N-2, 1); % 2nd order derivative of y_hat, D^2 y_hat

for i = 2: N-1
    D2yhat(i-1) = (y_hat(i+1) - 2*y_hat(i) + y_hat(i-1)) ./ dt^2;
end

SSE = SSE + lambda * (D2yhat' * D2yhat) * dt; 



