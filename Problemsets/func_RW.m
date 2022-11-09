function [Y] = func_RW(nbRW, T, sigma)

Y=nan(T,nbRW);
Y(1,:)=zeros(1,nb_process);

for j=1 : nbRW
    for t=2 : T
        Y(t,j) = Y(t-1, j) + randn() * sigma;
    end
end
