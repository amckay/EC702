DK = Grid.K/Kstar-1;
DKp = reshape(Kp,Grid.nK,Grid.nZ)./reshape(Grid.KK,Grid.nK,Grid.nZ) - 1;
plot(DK, DKp);
hold on;
plot(DK, zeros(Grid.nK,1), 'k--');
xlabel('K in % deviation from steady state')
ylabel('(K'' - K)/K')