%DEMO_CONDITIONING  How many digits survive?   (Math 6610, Lecture 8)

Agood = [2 1; 1 3];            % kappa ~ 3
Abad  = [1 1; 1 1+1e-8];       % kappa ~ 4e8

for A = {Agood, Abad}
    A = A{1};
    y = [1; 1];  b = A*y;                    % exact data, exact answer
    d = 1e-14 * norm(b) * [1; -1]/sqrt(2);   % tiny perturbation of b
    yp = A\(b + d);

    fprintf('kappa = %8.1e | input err = %.1e | output err = %.1e\n', ...
            cond(A), norm(d)/norm(b), norm(yp-y)/norm(y));
end
% kappa =  2.6e+00 | input err = 1.0e-14 | output err = 2.5e-14  (nothing lost)
% kappa =  4.0e+08 | input err = 1.0e-14 | output err = 4.0e-06  (8 digits gone)
