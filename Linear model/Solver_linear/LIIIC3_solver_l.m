function [t, u] = LIIIC3_solver_l(A, f, u0, tspan, h)

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

Ark = [1/6, -1/3, 1/6;
       1/6, 5/12, -1/12;
       1/6, 2/3, 1/6];
b = [1/6; 2/3; 1/6];
c = [0; 0.5; 1];

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);

    t1 = tn + c(1)*curr_h;
    t2 = tn + c(2)*curr_h;
    t3 = tn + c(3)*curr_h;

    f1 = f(yn, t1);
    f2 = f(yn, t2);
    f3 = f(yn, t3);

    M = sparse(3*N, 3*N);

    M(1:N, 1:N) = speye(N) + curr_h * Ark(1,1) * A;
    M(1:N, N+1:2*N) = curr_h * Ark(1,2) * A;
    M(1:N, 2*N+1:end) = curr_h * Ark(1,3) * A;

    M(N+1:2*N, 1:N) = curr_h * Ark(2,1) * A;
    M(N+1:2*N, N+1:2*N) = speye(N) + curr_h * Ark(2,2) * A;
    M(N+1:2*N, 2*N+1:end) = curr_h * Ark(2,3) * A;

    M(2*N+1:end, 1:N) = curr_h * Ark(3,1) * A;
    M(2*N+1:end, N+1:2*N) = curr_h * Ark(3,2) * A;
    M(2*N+1:end, 2*N+1:end) = speye(N) + curr_h * Ark(3,3) * A;

    rhs = repmat(yn, 3, 1);
    rhs(1:N) = rhs(1:N) + curr_h * (Ark(1,1)*f1 + Ark(1,2)*f2 + Ark(1,3)*f3);
    rhs(N+1:2*N) = rhs(N+1:2*N) + curr_h * (Ark(2,1)*f1 + Ark(2,2)*f2 + Ark(2,3)*f3);
    rhs(2*N+1:end) = rhs(2*N+1:end) + curr_h * (Ark(3,1)*f1 + Ark(3,2)*f2 + Ark(3,3)*f3);

    U = M \ rhs;
    U1 = U(1:N);
    U2 = U(N+1:2*N);
    U3 = U(2*N+1:end);

    K1 = -A * U1 + f1;
    K2 = -A * U2 + f2;
    K3 = -A * U3 + f3;

    u(:,n+1) = yn + curr_h * (b(1)*K1 + b(2)*K2 + b(3)*K3);
end

end
