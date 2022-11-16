function [nlag] = LagOrderSelectionARp(y,const,pmax,crit)

val_crit = nan(pmax,1);

% n : nb de parametres estimés : const + nb de phi
for lag = 1 : pmax
    Teff = size(y,1) - pmax;
    n = lag + const;
    sigutilde = ARpML(y,lag,const,0.05).sigutilde;
    if crit == "AIC"
        val_crit(lag) = log(sigutilde^2) + (2/Teff) * n;
    elseif crit == "SIC"
        val_crit(lag) = log(sigutilde^2) + (log(Teff)/Teff) * n;
    elseif crit == "HQC"
        val_crit(lag) = log(sigutilde^2) + (2*log(log(Teff))/Teff) * n;
    end
end
clc;


nlag = find(val_crit == min(val_crit));

end