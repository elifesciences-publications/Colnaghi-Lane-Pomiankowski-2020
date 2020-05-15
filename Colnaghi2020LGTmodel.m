%%%. Simulates Muller's ratchet in an asexual population or in a population
%%%. that undergoes LGT at a rate l.
%%%.
%%%. Population size = N, number of genes = g; evolve for nGen generations.
%%%.
%%%. The state of the system is described by a N x g matrix whose elements
%%%. correspond to the number of mutations at a certain locus
%%%.
%%%. Strength of selection s, fitness = (1 - s)^k where k total number of
%%%. mutations.
%%%. Mutation rate/genome/generation: U
%%%.
%%%. This script calculates the number of mutations that reach fixation 
%%%. after nGen generations.

figure;hold on
colors=get(gca,'colororder');

%%%. System Parameters
nGen=5000;                   % Number of generations
N=10000;                     % Population size
g=1000;                      % Number of genes; 1 gene =~ 1000bp
tic
nRuns = 50;                  % # of Iterations
Lvector = [0 1 5 10 g/5];    % eDNA length
l = 0.01;                    % LGT rate / individual / generation
U=0.01;                      % Genome-wide mutation rate
s = 0.001;                   % Strength of seleciton agains deleterious mutations
inLoad = 0.0;                % Initial mutation load


%%%. Simulation
disp(['n_0 = ',num2str(1/exp(U/s) *N)])
dmdt=zeros(numel(Lvector),nRuns);

for L1 = 1:numel(Lvector)           
    L=Lvector(L1);
    n_fixedM=zeros(nRuns,nGen);                 % # Fixed mutations
    totFixMut=zeros(1,nRuns);
    
    for r1 = 1:nRuns
        %%%. Initialise
        rand_mat = rand(N,g);                 % generate random matrix (for mutations)
        M = double(rand_mat < inLoad);        % N x g matrix; 0 = wild-type, 1 = mutant

        %%%. Evolution
        for t=2:nGen
            oldMat=M;
            M = offspring(M,s);                 % next generation
            M = mutate(M,U);                    % introduce new  mutations
 
            if l>0
                M = LGT2(M,l,L,oldMat);         % LGT
            end
 
            totFixMut(r1) = sum(min(M));
        end
        
    disp(['run ',num2str(r1),' of ', num2str(nRuns), '. L = ', num2str(L), '. g = ', num2str(g), '. *** ', num2str(totFixMut(r1)), ' fixed mutations ***'])
    end
    
    dmdt(L1,:) = totFixMut;
    bplot(dmdt(L1,:)/nGen,L1,'color',colors(L1,:))
    drawnow
end
toc











%%%. Additional functions
function X = offspring(X,s)
 N = numel(X(:,1));                            % population size 
 nMut = sum(X,2);                              % individual mutation load
 sPower = (1-s).^nMut;
 mFitness = mean(sPower);
 fitVec = sPower/mFitness;       % relative fitness
 y= randsample(N,N,true,fitVec);  
 X = X(y,:);
end


%%% Mutation function. It draws a random integer X from a Poisson
%%% distribution for each individual, which is the number of new mutations
%%% acquired. For each individual, it selects X loci in its genome which,
%%% if unmutated, acquire a new mutation.

function X = mutate(X,U)
 N = numel(X(:,1));                            % population size 
 g = numel(X(1,:));                            % genome size 
 mutString = random('Poisson',U,N,1);

 nMut = sum(mutString);
 if nMut>0
 yStr = randi([1,g],1,nMut); a=1;
 pStr = zeros (1,nMut);
 mut=find(mutString>0);
 for k=1:numel(mut)
  i=mut(k);
  x=repmat(i,1,mutString(i));
  pStr(a:a+mutString(i)-1)=x;
  a=a+mutString(i);
 end
  z= sum ([g.*(pStr-1);yStr]);
  X(z) = X(z)+1;
 end
  
end


%%% LGT function. It takes as input the probability of lgt l.
%%% string length L and the old matrix.
%%%
%%% For each individual, generate a random number; if it is smaller than l,
%%% select a random locus x between 1 and g-L+1. Select a random row from the
%%% old matrix. Elements x to x+L-1 of the new matrix become equal to those
%%% of the previous matrix. 

function X = LGT2(X,l,L,oldMat)
 N = numel(X(:,1));                            % population size 
 g = numel(X(1,:));                            % genome size 
 
 nLGT = random('Bino',N,l);                    % # of LGT transfer events in the population
 whoLGT = randsample(1:N,nLGT);                 % who undergoes LGT

 if sum(nLGT)>0
    LGTloci = randi([1 g-L+1],1,sum(nLGT));   % locus where rec begins
    X(whoLGT, LGTloci:LGTloci+L-1) = oldMat(whoLGT, LGTloci:LGTloci+L-1);
                                                % recombination 
 end
end
