%DEMO_DAMPING  Two test cases for dampnewton.m   (Math 6610, Lecture 8)

%% Test 1 -- scalar problem where the full Newton step diverges
f  = @(x) atan(x);
Df = @(x) 1/(1+x^2);

[~, hd, ld] = dampnewton(f, Df, 2);        % damped: converges
xp = 2;                                     % plain Newton, for comparison
for k = 1:5, xp(k+1) = xp(k) - atan(xp(k))*(1+xp(k)^2); end

fprintf('plain  Newton: |x^5| = %.2e\n', abs(xp(end)));   % 2.34e+10
fprintf('damped Newton: |x^5| = %.2e\n', abs(hd(end)));   % 2.71e-18
fprintf('damping used : '); fprintf('%g ', ld); fprintf('\n');  % 0.5 1 1 1 1

figure
semilogy(0:numel(xp)-1, abs(xp), 's--', 0:size(hd,1)-1, abs(hd), 'o-')
xlabel('iteration k'); ylabel('|x^k|')
legend('Newton (\lambda = 1)', 'damped Newton', 'Location', 'southwest')

%% Test 2 -- circle and parabola (Lecture 5), from an awkward start
f2  = @(x) [x(1)^2 + x(2)^2 - 1;  x(2) - x(1)^2];
Df2 = @(x) [2*x(1), 2*x(2);  -2*x(1), 1];

[x2, h2, l2] = dampnewton(f2, Df2, [0.02; 4]);
fprintf('first damping factors: '); fprintf('%g ', l2(1:5)); fprintf('\n');
fprintf('largest |x| visited  : %.1f   (plain Newton reaches 47)\n', max(abs(h2(:))));
