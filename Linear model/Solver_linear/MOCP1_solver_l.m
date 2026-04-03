function [t, u] = MOCP1_solver_l(A, f, dfdu, u_init, tspan, h)

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

    f_new = f(u(:,n), t_new);

    M = alpha2 * speye(N) + curr_h * beta2 * A;
    b = curr_h * beta2 * f_new - past_contrib;
    U = M \ b;
    u(:,n+1) = U;

    history = [history(:,2), U];
    t_hist = [t_hist(2), t_new];
end

end
