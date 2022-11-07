function [r_k] = func_ACFPlots(y, pmax, alph)
T = size(y,1); % get number of observations
y_demeaned = y - mean(y); % y in deviation from mean
r_k = nan(1, pmax); % Initialize output vector

% Compute variance
c0 = (1/T)*(y_demeaned'*y_demeaned);

for k=1 : pmax
    ck = (1/T)*y_demeaned(1+k:T,:)'*y_demeaned(1:T-k,:);
    r_k = ck/c0;
end

critval = norminv(1-alph/2);
ul = repmat(critval/sqrt(T), pmax, 1);
ll = -1*ul;

figure('name', 'Autocorrelation');
bar(r_k);
hold on;
plot(1:pmax, ul, 'Color','black');
plot(1:pmax, ll, 'Color','black');
hold off;


end