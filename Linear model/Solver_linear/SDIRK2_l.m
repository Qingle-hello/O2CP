function [t, u] = SDIRK2_l(A, f, dfdu, u0, tspan, h)

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
I = speye(N);

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);

    % Stage 1
    Z = yn;
    t_stage = tn + gamma * curr_h;
    g1 = f(yn, t_stage);
    RHS = Z + curr_h * gamma * g1;
    M = I + curr_h * gamma * A;
    U1 = M \ RHS;
    K1 = -A * U1 + g1;

    % Stage 2
    Z = yn + curr_h * (1 - gamma) * K1;
    t_stage = tn + curr_h;
    g2 = f(yn, t_stage);
    RHS = Z + curr_h * gamma * g2;
    U2 = M \ RHS;
    K2 = -A * U2 + g2;

    % Update
    u(:,n+1) = yn + curr_h * ((1 - gamma) * K1 + gamma * K2);
end

end
