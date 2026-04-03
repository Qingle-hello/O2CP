function [t, u] = BDF6_solver(A, f, dfdu, u_init, tspan, h)

if nargin < 3 || isempty(dfdu)
    has_jac = false;
    dfdu = @(u,t) sparse(length(u), length(u));
else
    has_jac = true;
end

N = size(u_init,1);
if size(u_init,2) ~= 6
    error('u_init must be N x 6 matrix.');
end

t0 = tspan(1); tf = tspan(2);
t = t0:h:tf;
if t(end) ~= tf
    t(end) = tf;
    warning('Last step adjusted.');
end
nt = length(t);
u = zeros(N, nt);
u(:,1) = u_init(:,6);

alpha_past = [-120/49, 150/49, -400/147, 75/49, -24/49, 10/147];
beta = 20/49;
coeffs = alpha_past(end:-1:1)';

history = u_init;

for n = 1:nt-1
    curr_h = t(n+1) - t(n);
    sum_past = history * coeffs;
    U0 = u(:,n);
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 100, 'TolFun', 1e-15);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    t_new = t(n+1);
    [U, ~, exitflag] = fsolve(@(U) resfun(U, sum_past, curr_h, t_new, A, f, dfdu, beta, has_jac, N), U0, options);
    if exitflag <= 0
        warning('fsolve failed at step %d.', n);
    end
    u(:,n+1) = U;
    history = [history(:,2:6), U];
end

end

function [res, jac] = resfun(U, sum_past, h, t_new, A, f, dfdu, beta, has_jac, N)
K = -A * U + f(U, t_new);
res = U + sum_past - h * beta * K;
if nargout > 1 && has_jac
    Df = dfdu(U, t_new);
    jac = speye(N) + h * beta * (A - Df);
else
    jac = [];
end
end
