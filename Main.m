clear;
addpath(genpath(pwd));
load('FuncAnaData.mat')
 
% data
chi = af(:, 1);
xi = af(:, 2);
y = yt(:, 1);

% parameters
lambda1 = 2;
lambda2 = 0.02;
deltat = dt; 
L1 = 16; 
L2 = 10; 
m1 = 3;
m2 = 3; 
ifPlot = true;

% B-spline basis
N = length(chi);
n1 = ceil(N / L1);
if (L1-1)*n1 > N
    n1 = n1-1;
end
n2 = ceil(N / L2);
if (L2-1)*n2 > N
    n2 = n2-1;
end
knot1 = [repelem(1, m1), (1: L1-1) * n1, repelem(N, m1)];
knot2 = [repelem(1, m2), (1: L2-1) * n2, repelem(N, m2)];
B1 = zeros(N, L1+m1-1); % B_{i,m1}^{(1)}, the basis for functions chi(t) and xi(t)
B2 = zeros(N, L2+m2-1); % B_{j,m2}^{(2)}, the basis for functions alpha1(t) and alpha2(t)

for i = 1: L1+m1-1
    B1(:, i) = B_Spline(i, knot1, m1);
end

for i = 1: L2+m2-1
    B2(:, i) = B_Spline(i, knot2, m2);
end

% stage 1: functional analysis for chi and xi 
D2basis1 = zeros(N-2, L1+m1-1); % second order derivative of phi 
for i = 2: N-1
    D2basis1(i-1, :) = (B1(i+1, :) - 2*B1(i, :) + B1(i-1, :))./ deltat^2 ;
end
R1 = D2basis1' * D2basis1 * deltat; 

c_hat_chi = inv(B1' * B1 + lambda1 * R1) * B1' * chi;
c_hat_xi = inv(B1' * B1 + lambda1 * R1) * B1' * xi;

chi_hat = B1 * c_hat_chi; % chi(t)
xi_hat = B1 * c_hat_xi; % xi(t)
c_hat = [c_hat_chi, c_hat_xi];

if ifPlot
    figure;
    hold on;
    plot(1: N, chi, "k");
    plot(1: N, chi_hat, "r");
    legend(["\chi", "Estimated \chi"]);
    hold off;

    figure;
    hold on;
    plot(1: N, xi, "k");
    plot(1: N, xi_hat, "r");
    legend(["\xi", "Estimated \xi"]);
    hold off;
end

% stage 2: functional regression
% constraints
A = []; % A*x <= b
b = [];
Aeq = []; % Aeq*x = beq
beq = [];
lb = []; % lb <= x <= ub 
ub = [];
%c = @(x) []; % c(x) <= 0 
%ceq = @(x) []; % ceq(x) = 0
%nonlcon = @(x) deal(c(x), ceq(x));
nonlcon = @(x) Fix1End(x, y, chi, xi, lambda2, B1, c_hat, B2, deltat, 2);

func_sse = @(c) SSE_Penalty_FR(c, y, chi, xi, lambda2, B1, c_hat, B2, deltat);
par0 = repelem(1, 2*(L2+m2-1)+3);
options = optimoptions(@fmincon, 'TolFun', 1e-06, 'TolX', 1e-06, 'MaxIter',100000, 'MaxFunEvals', 200000);
[coe, fval, ~] = fmincon(func_sse, par0, A, b, Aeq, beq, lb, ub, nonlcon, options);

[SSE, y_hat] = SSE_Penalty_FR(coe, y, chi, xi, lambda2, B1, c_hat, B2, deltat);

if ifPlot
    figure;
    hold on;
    plot(1: N, y, "k");
    plot(1: N, y_hat, "r");
    legend(["y", "y hat"]);
    hold off;
end


