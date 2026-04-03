function [t, u] = BE(A, f, dfdu, u0, tspan, h)

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

Ark = 1;
b = 1;
c = 1;

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = yn;
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 100, 'TolFun', 1e-15);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    [U, ~, exitflag] = fsolve(@(U) resfun(U, yn, curr_h, tn, A, f, dfdu, Ark, has_jac, N, c), U0, options);
    if exitflag <= 0
        warning('fsolve failed at step %d.', n);
    end
    t1 = tn + c*curr_h;
    K1 = -A * U + f(U, t1);
    u(:,n+1) = yn + curr_h * (b*K1);
end

end

function [res, jac] = resfun(U, yn, h, tn, A, f, dfdu, Ark, has_jac, N, c)
t1 = tn + c*h;
K1 = -A * U + f(U, t1);
res = U - yn - h * (Ark*K1);
if nargout > 1 && has_jac
    Df1 = dfdu(U, t1);
    jac = speye(N) + h * Ark * (A - Df1);
else
    jac = [];
end
end
