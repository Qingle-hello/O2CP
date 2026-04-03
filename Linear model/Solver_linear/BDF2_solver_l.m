function [t, u] = BDF2_solver_l(A, f, dfdu, u_init, tspan, h)

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

alpha_past = [-4/3, 1/3];
beta = 2/3;
coeffs = alpha_past(end:-1:1)';

history = u_init;
I = speye(N);

for n = 1:nt-1
    curr_h = t(n+1) - t(n);
    sum_past = history * coeffs;
    t_new = t(n+1);
    g = f(u(:,n), t_new);
    RHS = -sum_past + curr_h * beta * g;
    M = I + curr_h * beta * A;
    U = M \ RHS;
    u(:,n+1) = U;
    history = [history(:,2), U];
end

end
