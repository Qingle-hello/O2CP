function [t, u] = OCP(A, f, dfdu, u0, tspan, h)

a0=1; a1=-0.21014; a2=0.00486; eb1=0.78986; eb2=0.38283;
p0=1; p1=0.37797;

coeff = [eb2, eb1, a0];
roots_lambda = roots(coeff);
c = coeff(1); 
x1 = -roots_lambda(1); 
x2 = -roots_lambda(2);

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

I = spdiags(ones(N,1),0,N,N);

for n = 1:nt-1
    yn = u(:,n);
    tn = t(n);
    curr_h = t(n+1) - t(n);
    U0 = yn;

    u(:,n+1) = real( ...
                a2/eb2 * U0 + ...
             (c*(x1*I + curr_h*A))\((x2*I + curr_h*A)*...
               (((a0*eb2-a2)/eb2*I + (a1*eb2-a2*eb1)/eb2 * curr_h*A) * U0 + ...
               curr_h * ((p0*I + p1 * curr_h*A)*f(U0,tn)))) ...
                   );
end

end
