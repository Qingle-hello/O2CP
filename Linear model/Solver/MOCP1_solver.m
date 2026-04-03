function [t, u] = MOCP1_solver(A, f, dfdu, u_init, tspan, h)

if nargin < 3 || isempty(dfdu)
    has_jac = false;
    dfdu = @(u,t) sparse(length(u), length(u));
else
    has_jac = true;
end

N = size(u_init,1);
if size(u_init,2) ~= 2
    error('u_init must be N x 2 matrix.');
end

t0 = tspan(1); tf = tspan(2);
t = t0:h:tf;
if t(end) ~= tf
    t(end) = tf;
    warning('Last step adjusted.');
end
nt = length(t);
u = zeros(N, nt);
u(:,1) = u_init(:,2);

alpha2 = 1;
beta2 = 0.56380;
alpha_hist = [-0.02178, -0.97822];
beta_hist = [0.00047, 0.46300];

history = u_init;
t_hist = [t0 - h, t0];

for n = 1:nt-1
    curr_h = t(n+1) - t(n);
    t_new = t(n+1);

    past_contrib = zeros(N,1);
    for i = 1:2
        ui = history(:,i);
        ti = t_hist(i);
        Ki = A * ui - f(ui, ti);
        past_contrib = past_contrib + alpha_hist(i) * ui + curr_h * beta_hist(i) * Ki;
    end

    U0 = u(:,n);
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 10, 'TolFun', 1e-16);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    [U, ~, exitflag] = fsolve(@(U) resfun(U, past_contrib, curr_h, t_new, A, f, dfdu, beta2, alpha2, has_jac, N), U0, options);
    if exitflag <= 0
        warning('fsolve failed at step %d.', n);
    end
    u(:,n+1) = U;

    history = [history(:,2), U];
    t_hist = [t_hist(2), t_new];
end

end

function [res, jac] = resfun(U, past_contrib, h, t_new, A, f, dfdu, beta2, alpha2, has_jac, N)
K = A * U - f(U, t_new);
res = alpha2 * U + h * beta2 * K + past_contrib;
if nargout > 1 && has_jac
    Df = dfdu(U, t_new);
    jac = alpha2 * speye(N) + h * beta2 * (A - Df);
else
    jac = [];
end
end
