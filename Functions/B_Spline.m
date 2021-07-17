function [B] = B_Spline(t, knot, m)

% Calculate the value of the t-th B-spline basis function
% Inputs: 
%   t: index of basis, 1 <= t <= m+L-1
%   knot: sequence of L+2m-1 knots
%   m: the order of polynomials, which is one more than its degree
% Outputs: 
%   B: value of B-spline basis function 

B = zeros(1, knot(end)); % B_{t,m}
time = 1: knot(end);

if m > 1
    B1 = B_Spline(t, knot, m-1); % B_{t,m-1}
    B2 = B_Spline(t+1, knot, m-1); % B_{t+1,m-1}
    if knot(t+m-1) - knot(t) ~= 0
        alpha1 = (time - knot(t)) ./ (knot(t+m-1) - knot(t)); % alpha_{t,m}
    else
        alpha1 = 0;
    end
    
    if knot(t+m) - knot(t+1) ~= 0
        alpha2 = 1 - (time - knot(t+1)) ./ (knot(t+m) - knot(t+1)); % 1 - alpha_{t+1,m}
    else
        alpha2 = 0;
    end
    
    B = alpha1 .* B1 + alpha2 .* B2; % B_{t,m} = alpha_{t,m} * B_{t,m-1} + [1 - alpha_{t+1,m}] * B_{t+1,m-1}
else
    if knot(t+1) ~= knot(end)
        B(knot(t): knot(t+1)-1) = 1;
    else
        B(knot(t): end) = 1;
    end
end









