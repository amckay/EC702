
T = 1000;
Sim = Simulate(Par,bKp,'random',T);


% Compute moments from simulated data
disp(['St. dev. of log Y = ' num2str(std(log(Sim.Y))) ]);
disp(['St. dev. of log C = ' num2str(std(log(Sim.C))) ]);
CY_corr = corrcoef(Sim.C,Sim.Y);
disp(['Corr. of C and Y = ' num2str(CY_corr(1,2)) ]);