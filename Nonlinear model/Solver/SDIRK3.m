function [t, u] = SDIRK3(A, f, dfdu, u0, tspan, h)
 
if nargin < 3 || isempty(dfdu)
    has_jac = false;
    dfdu = @(u,t) sparse(length(u), length(u));
else
    has_jac = true;
end
N = length(u0);
t0 = tspan(1); tf = tspan(2);
t = t0:h:tf;
if t(end) ~= tf
    t(end) = tf;
    warning('Last step adjusted.');
end
nt = length(t);
u = zeros(N, nt);
u(:,1) = u0;
gamma = 0.4358665215;
tau = (1 + gamma) / 2;
b1 = -(6 * gamma^2 - 16 * gamma + 1) / 4;
b2 = (6 * gamma^2 - 20 * gamma + 5) / 4;
c = [gamma, tau, 1];
a = [gamma, 0, 0; ...
     tau - gamma, gamma, 0; ...
     b1, b2, gamma];
b = [b1, b2, gamma];
for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = yn;
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 10, 'TolFun', 1e-12);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    K = zeros(N, 3);
    U = zeros(N, 3);
    for i = 1:3
        Z = yn;
        for j = 1:i-1
            Z = Z + curr_h * a(i,j) * K(:,j);
        end
        Ark = a(i,i);
        t_stage = tn + c(i) * curr_h;
        [U(:,i), ~, exitflag] = fsolve(@(Ui) resfun(Ui, Z, curr_h, t_stage, A, f, dfdu, Ark, has_jac, N), U0, options);
        if exitflag <= 0
            warning('fsolve failed at stage %d, step %d.', i, n);
        end
        K(:,i) = -A * U(:,i) + f(U(:,i), t_stage);
        U0 = U(:,i);   
    end
    u(:,n+1) = yn + curr_h * (b(1)*K(:,1) + b(2)*K(:,2) + b(3)*K(:,3));
end
end

function [res, jac] = resfun(U, Z, h, t_stage, A, f, dfdu, Ark, has_jac, N)
K = -A * U + f(U, t_stage);
res = U - Z - h * Ark * K;
if nargout > 1 && has_jac
    Df = dfdu(U, t_stage);
    jac = speye(N) + h * Ark * (A - Df);
else
    jac = [];
end
end
