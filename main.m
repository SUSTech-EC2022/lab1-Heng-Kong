fitFunc = 'objFunc'; % name of objective/fitness function
T = 200; % total number of evaluations
args.min = 0;
args.max = 31;
args.pop_size = 10;
args.gen_len = 5;
[bestSoFarFit, bestSoFarSolution] = simpleEA(fitFunc, T,args);