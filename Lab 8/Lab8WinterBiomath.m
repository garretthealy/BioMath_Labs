%% Lab 8: Gillespie simulations of continuous-time Markov processes
%
% Name: Garrett Healy
%
% This algorithm is used to simulate the behavior of a CTMC where there are
% different species of individuals and random transitions can be described
% by *reactions* that modify the numbers of individuals in at least one
% species. Before implementing a Gillespie simulation, write down all
% possible reactions and how they affect each species, and assign each one
% an index.
%
% Gillespie algorithm
%
% # Set values for the rates c(i) where i is the reaction index and is between 1 and M
% # Set initial states x(j), where j is the index of each species and is between 1 and N
% # Set time t to 0, and set maximum time Tmax
% # while (t<Tmax):
% # compute propensities a(i) for each reaction i (according to the rules of the reactions) and set A=sum(a)
% # generate uniform random numbers r1 and r2 between 0 and 1
% # compute the time of the next reaction tau=-log(r1)/A
% # the next reaction index k is the integer such that r2*A is *between* the sum of a(i) from 1 to k-1 and the sum of a(i) from 1 to k
% # update the numbers of species x(j) according to which reaction took place
% # update time t = t+tau
%
% Grading: 
% Part 1: 1.1: 4 pts (3 for code, 1 for plots), 1.2: 4 pts (3 for code, 1 for plots), 
% 1.3: 3 pts (1 for code, 1 for extinction time calculation, 1 for comparison with ODE solution)
% Part 2: 2.1: 3 pts (1 for null clines, 1 for fixed pts, 1 for prediction of dynamics), 
% 2.2: 5 pts (3 pts for code, 1 for phase plane plot, 1 for answer),
% 2.3: 4 pts (3 pts for code, 1 for answer),
% 2.4: 5 pts (3 pts for code, 1 for variance plot, 1 for answer)
% 2.5: 6 pts (3 pts for code, 1 for each plot, 1 for answer describing the period and amplitude of oscillations in the two simulations)
% Total: 34 pts
%% Part 1: Birth and death processes
% 1.1 *Pure birth process*
%
% A pure birth process is a population where each individual produces one
% offspring with a given rate, independent of the others. In the language
% of the master equation, there is only one reaction possible: increasing
% the size of the population by 1, with effective rate (propensity) N.
% Implement the Gillespie algorithm for a pure birth process, initially
% starting with population of 1, with birth rate of 0.1 per day, and run it
% for 100 days. (Hint: since there is only one reaction, you don't need to
% decide which reaction takes place, only the time of the next birth.)
% The mean of this process is governed by the ODE:
%    dN/dt = 0.1N, with N(0)=1
% Plot the solution of this ODE on the same plot as the 10 stochastic
% realizations of the pure birth process.

M = 1;
N = 1; 

c = 0.1; 
X = 1; 

t = 0; %sets initial time 
tplot = zeros(1); 
tmax = 100;
k = 2; 
xplot = ones(1); 
tau=1; 

while (t<tmax) && tau~=0
    r1 = rand;
    r2 = rand; 
    a = c*X;  
    A = sum(a); 
    tau = -log(r1)/A;
    t = t+tau; 
    X = X+a; 
    tplot(k) = t; 
    xplot(k) = X; 
    k = k+1; 
end

plot(tplot,xplot)

% 1.2 *Pure death process*
%
% Modify your code to model a pure death process with a stochastic death
% rate for each individual, starting with a population of 1000 individuals,
% with a death rate of 0.05 per day, and run 10 simulations for 100 days.
% Compute the solution for the mean number of individuals from the ODE:
%
%   dN/dt = -0.05N, with N(0)=1000
%
% Plot the solution of this ODE on the same plot as the 10 stochastic
% realizations of the pure birth process.
%
% 1.3 *Extinction time statistics*
% Run 100 simulations until all the individuals die and compute the
% mean time to extinction and the standard deviation. How well is it
% predicted by the deterministic mean solution (exponential decay) dipping
% below 1 individual?  
%
%% Part 2. Predator-prey stochastic model
% Now let us model a system with more than one reaction, namely a
% predator-prey system. Let us call the number of prey X and the number of
% predators Y . The three different processes are: 
%
% # Birth of one prey with rate c1 = 20 and propensity a1 = c1*X
% # Death of one predator with rate c2 = 10 and propensity a2 = c2*Y
% # Predation (death of one prey and birth of one predator) with rate c3
% =0.01 and propensity  a3 = c3*X*Y 
%
%
% 2.1 *ODE model*
%
% Write down the two-variable differential equation (Lotka-Volterra) based
% on the model above. Find the nullclines and the fixed points and predict
% the qualitative behavior of the system. This is the model for the mean
% number of predators and prey.
%
% 2.2 *Stochastic model*
%
% Implement the Gillespie algorithm for the three-reaction predator-prey
% model. Using the parameters above and initial conditions X0 = 500, Y0 =
% 500, propagate the system for a large number of steps (e.g. 100000)
% and plot the trajectory in phase space (plot(X,Y)). Describe whether it
% agrees with the qualitative behavior of the continuous ODE model, and how
% the stochasticity manifests itself.
%
% 2.3 *Extinction of species*
%
% Change the initial populations to be 10 each (of predatory and prey) so
% that you see extinction of one of the species in some simulations. Report
% the fraction of 200 simulations which go extinct in 1000 steps. Change
% the initial populations to 20 each, run 200 simulations, and again report
% the fraction that goes extinct in 1000 steps.
%
% 2.4 *Fluctuations around the equilibrium*
%
% Set your initial populations to the non-zero fixed point that you found in 2.1
% and re-run the simulation for 10000 steps. Because this is a stochastic
% process, we do not expect that the populations will remain at the fixed
% points as the ODE predicts. Using the results of your simulation, calculate the
% variance in the Euclidean distance from the fixed point across 10
% simulations as a function of time. Plot variance vs. t. What is the
% relationship? For a pure diffusion process, the variance will increase
% linearly as a function of time. Is that the case here? 
%
% 2.5 *Comparison of amplitude and frequencies of oscillations*
%
% This exercise will require that you run the simulation twice: Once with
% predator and prey initial populations equal to 500, and once with initial
% populations set at the coexistence fixed point, and run the simulation
% for 100000 steps. For each simulation, plot the predator population and
% the prey population as a function of time on the same plot. Using your
% plots, compare the amplitudes and frequencies of population oscillation
% when the initial populations are equal to 500 and when the initial
% populations are equal to the fixed points. How are they similar or
% different? Are amplitude and frequency independent of one another? (this
% is a hallmark of linear oscillations)

