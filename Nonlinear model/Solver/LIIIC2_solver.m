function [t, u] = LIIIC2_solver(A, f, dfdu, u0, tspan, h)

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

Ark = [1/2, -1/2;
       1/2, 1/2];
b = [1/2; 1/2];
c = [0; 1];

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = repmat(yn, 2, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 100, 'TolFun', 1e-16);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    [U, ~, exitflag] = fsolve(@(U) resfun(U, yn, curr_h, tn, A, f, dfdu, Ark, has_jac, N, c), U0, options);
    if exitflag <= 0
        warning('fsolve failed at step %d.', n);
    end
    U1 = U(1:N); U2 = U(N+1:end);
    t1 = tn + c(1)*curr_h; t2 = tn + c(2)*curr_h;
    K1 = -A * U1 + f(U1, t1);
    K2 = -A * U2 + f(U2, t2);
    u(:,n+1) = yn + curr_h * (b(1)*K1 + b(2)*K2);
end

end

function [res, jac] = resfun(U, yn, h, tn, A, f, dfdu, Ark, has_jac, N, c)
U1 = U(1:N); U2 = U(N+1:end);
t1 = tn + c(1)*h; t2 = tn + c(2)*h;
K1 = -A * U1 + f(U1, t1);
K2 = -A * U2 + f(U2, t2);
res1 = U1 - yn - h * (Ark(1,1)*K1 + Ark(1,2)*K2);
res2 = U2 - yn - h * (Ark(2,1)*K1 + Ark(2,2)*K2);
res = [res1; res2];
if nargout > 1 && has_jac
    Df1 = dfdu(U1, t1);
    Df2 = dfdu(U2, t2);
    jac = sparse(2*N, 2*N);
    jac(1:N,1:N) = speye(N) + h * Ark(1,1) * (A - Df1);
    jac(1:N,N+1:end) = h * Ark(1,2) * (A - Df2);
    jac(N+1:end,1:N) = h * Ark(2,1) * (A - Df1);
    jac(N+1:end,N+1:end) = speye(N) + h * Ark(2,2) * (A - Df2);
else
    jac = [];
end
end
