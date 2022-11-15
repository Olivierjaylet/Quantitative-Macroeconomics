function val_crit = LagOrderSelectionARp(y,const,pmax,crit)

val_crit = nan(pmax,1);

% n : nb de parametres estimés : const + nb de phi
for lag = 1 : pmax
    Teff = size(y,1) - pmax;
    n = lag + const;
    siguhat = ARpOLS(y,lag,const,0.05).siguhat;
    if crit == "AIC"
        val_crit(lag) = log(siguhat^2) + (2/Teff) * n;
    elseif crit == "SIC"
        val_crit(lag) = log(siguhat^2) + (log(Teff)/Teff) * n;
    elseif crit == "HQC"
        val_crit(lag) = log(siguhat^2) + (2*log(log(Teff))/Teff) * n;
    end
end


nlag = find(min(val_crit));

end