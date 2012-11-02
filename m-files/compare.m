A1 = load('../results/compare-1.txt');
A2 = load('../results/compare-2.txt');

A = A2(2:end, :); % the first value is biased

idxA1 = find(A1(:,1) > max(A2(:,1)));

A = [A ; A1(idxA1,:)];

n = A(:,1);
A = A(:,2:end);

[nu nui] = unique(n);

B = [];
nui = [0; nui]
for i=1:length(nu)

	idx = (nui(i)+1):nui(i+1);
	B = [B; sum(A(idx,:), 1)/length(idx)];
end

figure(1);
loglog(nu, B(:,1), 'r*-');
hold on;
loglog(nu, B(:,2), 'gs-');
loglog(nu, B(:,3), 'bd-');
loglog(nu, B(:,4), 'ko-');
loglog(nu, B(:,5), 'y.--');
loglog(nu, B(:,6), 'r+--');
loglog(nu, B(:,7), 'yo--');
loglog(nu, B(:,8), 'bs--');
loglog(nu, B(:,9), 'c^--');
grid on;
hold off;
xlabel('n');
ylabel('T [s]');
legend('MUL-fft', 'MUL-lma', 'DIV-inv', 'DIV-lda', 'INV', 'SQRT-inv', 'SQRT-div', 'SQRT_4-inv', 'SQRT_4-div');
