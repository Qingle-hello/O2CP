function [t, u] = BE_IMEX(A, f, dfdu, u0, tspan, h)

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
I = spdiags(ones(N,1),0,N,N);

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = yn;
    u(:,n+1) = (I + curr_h*A) \ (U0 + curr_h * f(U0, tn));
end

end
