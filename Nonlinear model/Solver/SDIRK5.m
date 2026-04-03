function [t, u] = SDIRK5(A, f, dfdu, u0, tspan, h)
 


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
gamma = 1/4;
c = [1/4, 3/4, 11/20, 1/2, 1];
a = [1/4, 0, 0, 0, 0; ...
     1/2, 1/4, 0, 0, 0; ...
     17/50, -1/25, 1/4, 0, 0; ...
     371/1360, -137/2720, 15/544, 1/4, 0; ...
     25/24, -49/48, 125/16, -85/12, 1/4];
b = [25/24, -49/48, 125/16, -85/12, 1/4];
for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = yn;
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 10, 'TolFun', 1e-12);
    if has_jac
        options = optimoptions(options, 'SpecifyObjectiveGradient', true);
    end
    K = zeros(N, 5);
    U = zeros(N, 5);
    for i = 1:5
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
        U0 = U(:,i);  % Use previous U as initial guess for next stage
    end
    u(:,n+1) = yn + curr_h * (b(1)*K(:,1) + b(2)*K(:,2) + b(3)*K(:,3) + b(4)*K(:,4) + b(5)*K(:,5));
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
