function [t, u] = SDIRK2(A, f, dfdu, u0, tspan, h)
 
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
gamma = (2 - sqrt(2))/2;
for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = yn;
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 10, 'TolFun', 1e-10);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    % Stage 1
    Z = yn;
    Ark = gamma;
    t_stage = tn + gamma * curr_h;
    [U1, ~, exitflag] = fsolve(@(U) resfun(U, Z, curr_h, t_stage, A, f, dfdu, Ark, has_jac, N), U0, options);
    if exitflag <= 0
        warning('fsolve failed at stage 1, step %d.', n);
    end
    K1 = -A * U1 + f(U1, t_stage);
    % Stage 2
    Z = yn + curr_h * (1 - gamma) * K1;
    Ark = gamma;
    t_stage = tn + curr_h;
    U0 = U1;
    [U2, ~, exitflag] = fsolve(@(U) resfun(U, Z, curr_h, t_stage, A, f, dfdu, Ark, has_jac, N), U0, options);
    if exitflag <= 0
        warning('fsolve failed at stage 2, step %d.', n);
    end
    K2 = -A * U2 + f(U2, t_stage);
    % Update
    u(:,n+1) = yn + curr_h * ((1 - gamma) * K1 + gamma * K2);
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
