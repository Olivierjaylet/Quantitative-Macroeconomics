%%
ENDO=importdata('C:\Users\User\OneDrive\Bureau\Quantitative-Macroeconomics\VAR\USOil.csv');

data = ENDO.data;
const = 1;
[T, K]  = size(data);     % sample size

pmax = 12;
nobs_eff = T - pmax;

INFO_CRIT = nan(pmax,1);  % initialize
pmin=1;
crit='AIC';


Y = transpose(data((pmax+1):T,:));
Zmax = transpose(lagmatrix(data,[1:pmax]));


if const == 1
    Zmax = [ones(1,T); Zmax];
elseif const == 2
    Zmax=[ones(1,T); 1:T; Zmax];
end


for p=pmin : pmax
    n=const+p*K;
    Z = Zmax(1:n,pmax+1:end);
    A = (Y*Z')/(Z*Z'); % OLS and Gaussian ML estimate
    U = Y-A*Z; % OLS and ML residuals
    UUt = U*U'; % sum of squared residuals
    SIGMLu = (UUt./nobs_eff); % ML: not adjusted for # of estimated coefficients
    phi = p * K^2 + K;
    if strcmp(crit,'AIC') % Akaike
        INFO_CRIT(p,:) = log(det(SIGMLu)) + (2/T)*phi;
    elseif strcmp(crit,'SIC') % Schwartz
        INFO_CRIT(p,:) = log(det(SIGMLu)) + (log(T)/T)*phi;
    elseif strcmp(crit,'HQC') % Hannan-Quinn
        INFO_CRIT(p,:) = log(det(SIGMLu)) + (2*log(log(T))/T)*phi;
    end
end
results = [transpose(1:pmax) INFO_CRIT];
nlag = find(INFO_CRIT == min(INFO_CRIT));