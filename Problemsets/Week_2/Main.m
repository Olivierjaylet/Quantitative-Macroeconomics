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
AR1_mat = func_AR1(phi, T, sigma);
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