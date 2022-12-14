%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Week 2 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DB.Nomics : list of data providers
% Our world in Data : research ideas and the data source
% Try official sources first, then large organisations
% Look at papers using data
%% Prepare datas
import_data;
Date = NorwayGDP.DATE;
GDP = NorwayGDP.DATA;

%%%%%%%% Ex 1) Visualizing Time Series Data %%%%%%%%
%% Plot quaterly growth rate
for j=1:2
    if j==1
        str_sample = '1978-Q1 to 2022-Q2';
        sample_start = 1;
        sample_end = length(GDP);
    elseif j==2
        str_sample = '1980-Q1 to 2012-Q4';
        sample_start = find(Date=='1980-Q1');
        sample_end = find(Date == '2012-Q4');
    end
    gdp_growth = (GDP((sample_start+1):sample_end) - GDP((sample_start):(sample_end-1))) ./ GDP((sample_start):(sample_end-1));
    
    dates_subsample = Date((sample_start+1):sample_end);
    figure('name', ['Plot of Q-on-Q growth in Norway', str_sample]);
    plot(dates_subsample, gdp_growth)
end


%% Histograms + fitted Normal distribution
figure;
tiledlayout(2,2);
nexttile();
histogram(gdp_growth, 10);
nexttile();
histfit(gdp_growth, 10);
nexttile();
histogram(gdp_growth, 30);
nexttile();
histfit(gdp_growth, 30);

%% Q-Q Plot
figure;
qqplot(gdp_growth);
title("Normal probability plot");

%% Boxplot 
figure;
boxplot(gdp_growth);
title("Boxplot");

%% Estimate of the average Growth rate and its Std
mean(GDP, 1);
std(GDP);

%%
%%%%%%%% Ex 3) Fundamental concepts of Univariate TSA %%%%%%%%

%% Plot 200 Observations
T=200;
eps = randn(T, 1);
sigma = 1;

Serie1 = nan(T,1);
for t=1:T
    Serie1(t,1) = eps(t) * sigma;
end

Serie2 = nan(T,1);
for t=3:T-2
    Serie2(t,1) = (1/5)*(eps(t-2) + eps(t-1) + eps(t) + eps(t+1) + eps(t+2));
end

figure; 
plot(Serie1);

figure;
plot(Serie2);

%% Plot an AR1 process
phi = [-0.8 0.7 0.9 1.01];
sigma = 1;
T = 200;
AR1_mat = func_AR1_1(phi, T, sigma);
%%

str_phi = ["\phi=-0.8", "\phi=0.7", "\phi=0.9", "\phi=1.01"];

figure('Name','AR Plots');
tiledlayout(length(phi),1);
for i=1 : length(phi)
    subplot(2,2,i);
    plot(AR1_mat(:,i));
    title(str_phi(i));
end

%% Random Walks
sigma = 1;
nbRW = 16;
T = 200;
RW_mat = func_RW(nbRW, T, sigma);

figure('name', 'Random Walks : y_t = y_t{-1} + /varepsilon');
sgtitle('Random walk');
for i=1 : nbRW
    subplot(4,4,i);
    plot(RW_mat(:,i));
end

%% Autocorrelation Function (ACF)
r_k = func_ACFPlots(AR1_mat, 3, 95);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Week 3 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
sig_eps = 0.5;
Y(1,:) = repmat(mu,1,B);
for b = 1:B
    epsi = sig_eps * randn(T,1);
    for t=2:T
     Y(t,b)=c + phi*Y(t-1,b) + epsi(t);
    end
end

mu_hat = mean(Y);
%var_Y = sig_eps^2/(1-phi^2); % analytical variance of AR(1)
var_Y = sig_eps^2/(1-phi)^2; % Correct standardized variance


% Standardization
Z = sqrt(T).*(mu_hat - mu)./sqrt(var_Y);

x = -5:0.1:5;
figure('name', 'Centrale Limit Theorems');
histogram(Z, 'Normalization','pdf');
hold on;
plot(x, normpdf(x), 'LineWidth',2);
title('Dependant Data (correct)');
ylim([0 0.45]);
hold off;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Week 4 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=func_AR1_2(0.95, 1000, 1, 1);

[thetas, SE_thetas, t_stats] = ARpOLS(Y,1,1);
figure;
scatter(y_hat, y, 'blue');
a= lsline;
set(a(1),'color','r');

%% Estimate AR(4) model from the AR4.csv file

csv = readtable('AR4.csv');
Y_AR4 = table2array(csv);

%%
OLS = ARpOLS(Y_AR4,1,1,0.05);

% After data manipulation we also can use the function fitlm(Y,y)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Week 4 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_AR1 = func_AR1_2(0.95, 1000, 1, 1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Information Criteria for AR(p) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

AIC= LagOrderSelectionARp(Y_AR1,1,10,"AIC");
%%
SIC = LagOrderSelectionARp(Y_AR4,1,10,"SIC");
%%
HQC = LagOrderSelectionARp(Y_AR4,1,10,"HQC");


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Portmanteau test for residual autocorrelation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data
csv = readtable('gnpdeflator.csv');
gnpdeflator = table2array(csv);
inflation = zeros(size(gnpdeflator,1)-1, 1);
for i =2 : size(gnpdeflator,1)
    inflation(i) = log(gnpdeflator(i,3)) - log(gnpdeflator(i-1,3));
end

pmax = 12;      % set max number of lags
const = 0;      % no constant
alph = 0.05;    % significance level

%% Compute AIC without time trend
phat_AIC = LagOrderSelectionARp(inflation,const,pmax,"AIC"); % with constant

%% Estimate OLS with and AR(phat)

OLS_phat = ARpOLS(inflation,phat_AIC,1, 0.05);

%% Estimate ML with and AR(1)

ML_AR1 = ARpML(inflation,1,1, 0.05);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap Confidence Interval for AR(1) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




