ENDO=importdata('C:\Users\User\OneDrive\Bureau\Quantitative-Macroeconomics\VAR\USOil.csv');
%%
const =1;
pmax = 12;
pmin = 1;
crit = "AIC";

[nlag, Z_mat, Y_test] = LagOrderSelectionVARp(ENDO.data,const,pmin,pmax,crit);

opt.const = 1;
VAR = VARReducedForm(ENDO.data,2,opt);

% %%
% data = ENDO.data;
% [T, K]  = size(data);
% p = 3;
% pmax = 12;
% nobs_eff = T - p; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = transpose(data((p+1):T,:));
% %%
% Z = transpose(lagmatrix(data,[1:p]));
% %%
% ZZ3 = Z(:,p+1:T); %remove initial observations
% % add deterministic terms if any
% if opt.const == 1
%     ZZ3= [ones(1,nobs_eff); ZZ3];
% elseif opt.const == 2
%     Z1=[ones(1,nobs_eff); (p):nobs_eff+(p-1); Z1];
% end
%%
data = ENDO.data;
[T, K]  = size(data);     % sample size
pmax = 12;

const=1;
INFO_CRIT = nan(pmax,1);  % initialize
crit='AIC';
for p=1 : pmax
    nobs_eff = T - p;
    Y = transpose(data((p+1):T,:));
    Z = transpose(lagmatrix(data,[1:p]));
    Z = Z(:,p+1:T);

    if const == 1
        Z= [ones(1,nobs_eff); Z];
    elseif const == 2
        Z=[ones(1,nobs_eff); (p):nobs_eff+(p-1); Z];
    end
    A= (Y*Z')/(Z*Z'); % OLS and Gaussian ML estimate
    U = Y-A*Z; % OLS and ML residuals
    UUt = U*U';   % sum of squared residuals
    SIGMLu = UUt/nobs_eff; % ML: not adjusted for # of estimated coefficients

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