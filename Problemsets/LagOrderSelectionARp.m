function nlag = LagOrderSelectionARp(y,const,pmax,crit)

val_crit = nan(pmax,1); % Initialize 

% n : number of estimated parameters : const + nb de phi
for lag = 1 : pmax
    Teff = size(y,1) - pmax;                            % effective sample size
    n = lag + const;                                    % number of freely estimated parameters
    resid = ARpOLS(y,lag,const,0.05).resid;             % OLS residual
    sigma2u = resid' * resid /(Teff-n);                 % ML estimate of OLS residual 
    if strcmp(crit, 'AIC')
        val_crit(lag) = log(sigma2u) + (2/Teff) * n;
    elseif strcmp(crit,'SIC')
        val_crit(lag) = log(sigma2u) + (log(Teff)/Teff) * n;
    elseif strcmp(crit, 'HQC')
        val_crit(lag) = log(sigma2u) + (2*log(log(Teff))/Teff) * n;
    end
end
clc;

results = [transpose(1:pmax) val_crit];
nlag = find(val_crit == min(val_crit));

 
[~,idx] = sort(results(:,2));
results = results(idx,:);
disp(results);
% Display summary of results
fprintf('*************************************************************\n');
fprintf('*** OPTIMAL ENDOGENOUS LAGS FROM %s INFORMATION CRITERIA ***\n',crit);
fprintf('*************************************************************\n'); 
disp(array2table(results,'VariableNames',['Lag',crit]))
% fprintf(' Optimal number of lags (searched up to %d lags):\n',pmax);
% fprintf(' %s Info Criterion:   %d\n',crit,nlag);
% fprintf('\n');


end