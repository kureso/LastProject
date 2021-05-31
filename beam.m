f = @(x) 75*9.81*exp(-x.^2/(2*0.01)); alpha = 0; beta = 0; gama = 0 ;delta = 0; n = 18; E = 9.4*10^6; Izz = 3;l = 10;X = [-l/2,l/2];
L = E*Izz*diffmat([n n+4],4,X);
vT = diffrow(n+4,0,-l/2,X); wT = diffrow(n+4,0,l/2,X); uT = diffrow(n+4,1,-l/2,X); sT = diffrow(n+4,1,l/2,X);
A = [L; vT; wT; uT; sT];
rhs = [gridsample(f,n); alpha; beta; gama; delta];
u = A\rhs;
tiledlayout(2,1)
nexttile
plot(chebfun(-u),'.-');ylim([-1 1]);
nexttile
plot(chebfun(-u),'.-');ylim([-0.001 0.001])
