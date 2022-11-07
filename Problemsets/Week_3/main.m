%%%%%% Law of Large Numbers %%%%%%
T = 10000;
u = zeros(T, 1);

a = -50; % lower_bound
b = 50; % upper_bound

uni_dist = unifrnd(a, b, T, 1);
norm_dist = normrnd(0, 1, T, 1);

list_dist = [uni_dist, norm_dist];

for i =1 : size(list_dist, 2)
    dist = list_dist(:,i);
for t = 1 : T
     u(t) = sum(dist(1:t)) / t;
end
    figure;
    plot(u);
end
%% LLN for AR1

sigma=1;
epsi = rand(T, 1)*sigma;
y = zeros(T, 1);
mean_y = zeros(T, 1);
phi = 0.9;

for t=2 : T+1
    y(t) = phi * y(t-1) + epsi(t-1);
    mean_y(t-1) = sum(y(1:t))/t;
end

figure; 
plot(mean_y);



%% Central Limit Theorem
B = 5000;
T = 10000;
phi = 0.8;
c = 3;
mu = c/(1-phi);
sig_eps = 0.4;
Y = func_AR1(phi, T, B, sig_eps);

mu_hat = mean(Y);
var_Y = sig_eps^2/(1-phi^2); % analytical variance of AR(1)
%var_Y = sig_eps^2/(1-phi)^2; % Correct standardized variance


% Standardization
Z = sqrt(T).*(mu_hat - mu)./sqrt(var_Y);

figure;
histogram(Z);
