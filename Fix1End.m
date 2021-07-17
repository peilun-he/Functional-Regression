function [c, ceq] = Fix1End(x, y, chi, xi, lambda2, B1, c_hat, B2, deltat, clamped)

% Non-linear constraints on the 1st and last points. 
% Inputs: 
%   x: a vector of parameters
%   y: data
%   chi: time series chi
%   xi: time series xi
%   lambda2: lambda2
%   B1: B_{i,m1}^{(1)}, the basis for functions chi(t) and xi(t)
%   c_hat: [c^{chi}, c^{xi}]
%   B2: B_{j,m2}^{(2)}, the basis for functions alpha1(t) and alpha2(t)
%   deltat: delta t
%   clamped: 0 -> no constraints
%            1 -> only fix function values
%            2 -> fix both function values and 1st order derivatives;
% Outputs: 
%   c: non-linear constraints c(x) <= 0
%   ceq: non-linear constraints ceq(x) = 0

if clamped == 0
    c = [];
    ceq = [];
elseif clamped == 1
    [~, y_hat] = SSE_Penalty_FR(x, y, chi, xi, lambda2, B1, c_hat, B2, deltat);
    c = [];
    ceq = [y_hat(1) - y(1), y_hat(end) - y(end)];  
elseif clamped == 2
    beta1 = x(2);
    beta2 = x(3);
    n = (length(x) - 3) / 2; % the length of c^{alpha1} and c^{alpha2}
    coe_func_chi = x(4: 3+n)'; % c^{alpha1}
    coe_func_xi = x(4+n: end)'; % c^{alpha2}
    c_hat_chi = c_hat(:, 1); % c^{chi}
    c_hat_xi = c_hat(:, 2); % c^{xi}

    derivative_int = (B1*c_hat_chi) .* (B2*coe_func_chi) + (B1*c_hat_xi) .* (B2*coe_func_xi); % the derivatives of two integrals, alpha1(t)*chi(t) + alpha2(t)*xi(t)     
    derivative1 = beta1*(chi(2)-chi(1))/deltat + beta2*(xi(2)-xi(1))/deltat - derivative_int(1); % derivative at 1st point
    derivativen = beta1*(chi(end)-chi(end-1))/deltat + beta2*(xi(end)-xi(end-1))/deltat - derivative_int(end); % derivative at last point

    [~, y_hat] = SSE_Penalty_FR(x, y, chi, xi, lambda2, B1, c_hat, B2, deltat);

    c = [];
    ceq = [y_hat(1) - y(1), y_hat(end) - y(end), derivative1 - (y(2) - y(1))/deltat, derivativen - (y(end) - y(end-1))/deltat];      
end



