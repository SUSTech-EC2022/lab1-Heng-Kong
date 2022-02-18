function [bestSoFarFit ,bestSoFarSolution ...
    ]=simpleEA( ...  % name of your simple EA function
    fitFunc, ... % name of objective/fitness function
    T, ... % total number of evaluations
    input) % replace it by your input arguments

% Check the inputs
if isempty(fitFunc)
  warning(['Objective function not specified, ''' objFunc ''' used']);
  fitFunc = 'objFunc';
end
if ~ischar(fitFunc)
  error('Argument FITFUNC must be a string');
end
if isempty(T)
  warning(['Budget not specified. 1000000 used']);
  T = '1000000';
end
eval(sprintf('objective=@%s;',fitFunc));
% Initialise variables
nbGen = 0; % generation counter
nbEval = 0; % evaluation counter
bestSoFarFit = 0; % best-so-far fitness value
bestSoFarSolution = NaN; % best-so-far solution
%recorders
fitness_gen=[]; % record the best fitness so far
solution_gen=[];% record the best phenotype of each generation
fitness_pop=[];% record the best fitness in current population 
%% Below starting your code

% Initialise a population
%% TODO
pop = randi([0,31],4,1);  %生成4个个体，0-31之间
bin_pop = dec2bin(pop); %10进制转换为2进制

% Evaluate the initial population
%% TODOs
fitness = objective(pop); %计算种群中每个个体的适应度
fitness_pop = [fitness_pop, max(fitness)] %获取最优适应度
if bestSoFarFit < max(fitness)
    [bestSoFarFit,i] = max(fitness);
    bestSoFarSolution = pop(i);
end

nbGen = nbGen + 1; %迭代次数
nbEval = nbEval + 4;  %计算适应度次数

fitness_gen = [fitness_gen, bestSoFarFit];
solution_gen = horzcat(solution_gen, bestSoFarSolution);

% Start the loop
while (nbEval<T)
    crossoverProb = fitness./sum(fitness); % roulette-wheel selection
    offspringGenes = [];
    for i = 1:2
        parentIndexes = [];
        for j = 1:2
            r = rand();
            for index = 1:4
                if r>sum(crossoverProb(1:index-1)) && r<=sum(crossoverProb(1:index))
                    break;
                end
            end
            parentIndexes = [parentIndexes, index];
        end
        crossoverPoint = randi(4);
        offspringGenes = [offspringGenes; [bin_pop(parentIndexes(1),1:crossoverPoint), bin_pop(parentIndexes(2),crossoverPoint+1:end)]];
        offspringGenes = [offspringGenes; [bin_pop(parentIndexes(2),1:crossoverPoint), bin_pop(parentIndexes(1),crossoverPoint+1:end)]];
    end
    mutationProb = 1/5;
    for i = 1:4
        isMutation = rand(1,5)<mutationProb;
        offspringGenes(i,isMutation) = dec2bin('1'-offspringGenes(i,isMutation))';
    end
    bin_pop = offspringGenes;
    population = bin2dec(bin_pop);
    fitness = objective(population);
    [A,index] = sort(fitness,1,'descend');
    fitness_pop=[fitness_pop,A(1)]
    for i = 1:4
        if fitness(i) > bestSoFarFit
            bestSoFarFit = fitness(i);
            bestSoFarSolution = population(i,:);
        end
    end
    nbEval = nbEval + 4;
    nbGen = nbGen + 1;
    fitness_gen=horzcat(fitness_gen, bestSoFarFit);
    solution_gen=horzcat(solution_gen, bestSoFarSolution);
end
bestSoFarFit
bestSoFarSolution

figure,plot(1:nbGen,fitness_gen,'b') 
title('Fitness\_Gen')

figure,plot(1:nbGen,solution_gen,'b') 
title('Solution\_Gen')

figure,plot(1:nbGen,fitness_pop,'b') 
title('Fitness\_Pop')
% Reproduction (selection, crossver)
%% TODO

% Mutation
%% TODO








