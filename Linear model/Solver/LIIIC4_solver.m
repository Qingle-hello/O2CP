function [t, u] = LIIIC4_solver(A, f, dfdu, u0, tspan, h)

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

sqrt5 = sqrt(5);
Ark = [1/12, -sqrt5/12, sqrt5/12, -1/12;
       1/12, 1/4, 1/6 - 7*sqrt5/60, sqrt5/60;
       1/12, 1/6 + 7*sqrt5/60, 1/4, -sqrt5/60;
       1/12, 5/12, 5/12, 1/12];
b = [1/12; 5/12; 5/12; 1/12];
c = [0; (5 - sqrt5)/10; (5 + sqrt5)/10; 1];

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = repmat(yn, 4, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 100, 'TolFun', 1e-15);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    [U, ~, exitflag] = fsolve(@(U) resfun(U, yn, curr_h, tn, A, f, dfdu, Ark, has_jac, N, c), U0, options);
    if exitflag <= 0
        warning('fsolve failed at step %d.', n);
    end
    U1 = U(1:N); U2 = U(N+1:2*N); U3 = U(2*N+1:3*N); U4 = U(3*N+1:end);
    t1 = tn + c(1)*curr_h; t2 = tn + c(2)*curr_h; t3 = tn + c(3)*curr_h; t4 = tn + c(4)*curr_h;
    K1 = -A * U1 + f(U1, t1);
    K2 = -A * U2 + f(U2, t2);
    K3 = -A * U3 + f(U3, t3);
    K4 = -A * U4 + f(U4, t4);
    u(:,n+1) = yn + curr_h * (b(1)*K1 + b(2)*K2 + b(3)*K3 + b(4)*K4);
end

end

function [res, jac] = resfun(U, yn, h, tn, A, f, dfdu, Ark, has_jac, N, c)
U1 = U(1:N); U2 = U(N+1:2*N); U3 = U(2*N+1:3*N); U4 = U(3*N+1:end);
t1 = tn + c(1)*h; t2 = tn + c(2)*h; t3 = tn + c(3)*h; t4 = tn + c(4)*h;
K1 = -A * U1 + f(U1, t1);
K2 = -A * U2 + f(U2, t2);
K3 = -A * U3 + f(U3, t3);
K4 = -A * U4 + f(U4, t4);
res1 = U1 - yn - h * (Ark(1,1)*K1 + Ark(1,2)*K2 + Ark(1,3)*K3 + Ark(1,4)*K4);
res2 = U2 - yn - h * (Ark(2,1)*K1 + Ark(2,2)*K2 + Ark(2,3)*K3 + Ark(2,4)*K4);
res3 = U3 - yn - h * (Ark(3,1)*K1 + Ark(3,2)*K2 + Ark(3,3)*K3 + Ark(3,4)*K4);
res4 = U4 - yn - h * (Ark(4,1)*K1 + Ark(4,2)*K2 + Ark(4,3)*K3 + Ark(4,4)*K4);
res = [res1; res2; res3; res4];
if nargout > 1 && has_jac
    Df1 = dfdu(U1, t1);
    Df2 = dfdu(U2, t2);
    Df3 = dfdu(U3, t3);
    Df4 = dfdu(U4, t4);
    jac = sparse(4*N, 4*N);
    jac(1:N,1:N) = speye(N) + h * Ark(1,1) * (A - Df1);
    jac(1:N,N+1:2*N) = h * Ark(1,2) * (A - Df2);
    jac(1:N,2*N+1:3*N) = h * Ark(1,3) * (A - Df3);
    jac(1:N,3*N+1:end) = h * Ark(1,4) * (A - Df4);
    jac(N+1:2*N,1:N) = h * Ark(2,1) * (A - Df1);
    jac(N+1:2*N,N+1:2*N) = speye(N) + h * Ark(2,2) * (A - Df2);
    jac(N+1:2*N,2*N+1:3*N) = h * Ark(2,3) * (A - Df3);
    jac(N+1:2*N,3*N+1:end) = h * Ark(2,4) * (A - Df4);
    jac(2*N+1:3*N,1:N) = h * Ark(3,1) * (A - Df1);
    jac(2*N+1:3*N,N+1:2*N) = h * Ark(3,2) * (A - Df2);
    jac(2*N+1:3*N,2*N+1:3*N) = speye(N) + h * Ark(3,3) * (A - Df3);
    jac(2*N+1:3*N,3*N+1:end) = h * Ark(3,4) * (A - Df4);
    jac(3*N+1:end,1:N) = h * Ark(4,1) * (A - Df1);
    jac(3*N+1:end,N+1:2*N) = h * Ark(4,2) * (A - Df2);
    jac(3*N+1:end,2*N+1:3*N) = h * Ark(4,3) * (A - Df3);
    jac(3*N+1:end,3*N+1:end) = speye(N) + h * Ark(4,4) * (A - Df4);
else
    jac = [];
end
end
