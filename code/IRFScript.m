
T = 42;
Sim = Simulate(Par,bKp,'irf',T);

% create a steady state structure with same fields as Sim
SS.K = Kstar;
SS.Z = 0;
SS.Y = f(Par,SS.K,SS.Z) - (1-Par.delta)*SS.K;
SS.C = SS.Y -Par.delta*SS.K;

Sim.Z = exp(Sim.Z);
SS.Z = exp(0);

figure;
T = length(Sim.Y);
flds = fieldnames(Sim);
for i = 1:numel(flds)
    subplot(2,2,i);
    plot(1:T, 100*(Sim.(flds{i})-SS.(flds{i}))/SS.(flds{i}) );
    title(flds{i});
end
