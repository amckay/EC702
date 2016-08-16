function Accuracy(Par,Grid,bC)
    TestGrid = Grid;
    TestGrid.nK = 200;

    TestGrid.K = linspace(Grid.K(1),Grid.K(end),TestGrid.nK);
    [Ztest,Ktest] =meshgrid(Grid.Z,TestGrid.K);
    TestGrid.KK = Ktest(:);
    TestGrid.ZZ = Ztest(:);


    C = PolyBasis(TestGrid.KK,TestGrid.ZZ) * bC;
    Kp = f(Par,TestGrid.KK,TestGrid.ZZ) - C;
    CEuler =  EulerRHS(Par,TestGrid,Kp,bC);

    plot(100*(TestGrid.K/Par.Kstar-1), reshape(  log10(abs ( CEuler./C-1 )), 200,7) )
    ylim([-7 -2])
    xlabel('K in % deviation from steady state')
    ylabel('Absolute Euler equation error, log base 10')