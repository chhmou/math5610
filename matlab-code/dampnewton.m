function [x, hist, lams] = dampnewton(f, Df, x0, tol, lmin)
%DAMPNEWTON  Damped Newton's method with the affine-invariant monotonicity test.
%
%   [x, hist, lams] = DAMPNEWTON(f, Df, x0, tol, lmin)
%
%   f, Df  function handles: residual f(x) and Jacobian Df(x)
%   x0     starting point
%   tol    stop when norm(dx) <= tol         (default 1e-12)
%   lmin   smallest damping factor allowed   (default 1e-4)
%
%   hist   one row per iterate,  lams  the accepted damping factors.
%
%   Math 6610, Lecture 8.

if nargin < 4 || isempty(tol),  tol  = 1e-12; end
if nargin < 5 || isempty(lmin), lmin = 1e-4;  end

x = x0(:);  lam = 1;  hist = x.';  lams = [];

for k = 1:100
    [L,U,P] = lu(Df(x));              % ONE factorization per Newton step
    dx = U \ (L \ (P*(-f(x))));       % Newton correction
    if norm(dx) <= tol, return, end   % stopping rule (Lecture 7)

    lam = min(1, 2*lam);              % try to grow back toward a full step
    while lam >= lmin
        xt   = x + lam*dx;                   % trial point
        dbar = U \ (L \ (P*(-f(xt))));       % simplified correction: reuse LU
        if norm(dbar) <= (1 - lam/2)*norm(dx)
            break                            % monotonicity test passed
        end
        lam = lam/2;                         % rejected: halve and retry
    end
    if lam < lmin
        error('dampnewton:stalled', ...
              'no acceptable step at k=%d: check Df, or the starting point', k);
    end

    x = x + lam*dx;
    hist(end+1,:) = x.';  lams(end+1) = lam;      %#ok<AGROW>
end
error('dampnewton:maxit', 'no convergence in 100 steps');
