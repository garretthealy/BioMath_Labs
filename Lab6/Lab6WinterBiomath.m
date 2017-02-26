%% Lab 6: Markov Chain Monte Carlo simulations
%
% Metropolis-Hastings Algorithm
%
% Points breakdown:
%
% 1.1 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 1.2 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 1.3 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 1.4 7 pts (3 pts for code, 1 pt for each plot, 1 pt for answer)
% 2.1 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 2.2 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 2.3 2 pts (1 pt for plot, 1 pt for answer)
% 3.1 5 pts (3 pts for code, 1 pt for plot, 1 pt for answer)
% 4.1 6 pts (4 pts for code, 1 pt for plot, 1 pt for answer)
% 4.2 6 pts (4 pts for code, 1 pt for plot, 1 pt for answer)
% Total: 56 points

% # Initialize the first state to X0 and set t = 0
% # Set the target probabilities P(X) (for all states X) and proposal
% probabilities q(i,j) of generating candidate state i from current state j.
% # Generate a candidate state Y using q(i,j) by drawing a random number u 
% from a uniform distribution on (0,1)
% # Draw another random number v; 
%   if v <= min(1,(P(Y)*q(X,Y))/(P(X)*q(Y,X))) then set X = Y 
%   else keep the same state
% # Set t=t+1 and repeat steps 3 through 5
%
%% Part 1: Single 2-state ion channel
% Simulate an ion channel switching between open and closed states, 
% with a number (e.g. 0) denoting the open state and another (e.g. 1) 
% denoting the closed state. Make the target fractions of open
% and closed states a changeable parameter, making sure they add up to 1. 
%
% 1.1 Use the simple proposal distribution: q(i,j) = 0.5 regardless of the states 
% i and j. Start with an open state. Run your simulation for 1000 time 
% steps using the target distribution: P(O) = 0.5; P(C) = 0.5. Plot the 
% states of the channel (0 and 1) over 1000 time steps. 
%
% 1.2 Repeat the simulation with the target distribution of P(O) = 0.1;
% P(C) = 0.9, and plot the states of the channel over 1000 time steps.
%
% 1.3 Change the proposal distribution so that q(O,C) = 0.99 and q(O,O) = 0.99 
% (i.e., the probability of proposing the open state is always 0.99) and 
% report if you see a difference in the simulation for both of the two
% target distributions.
%
% 1.4 Using P(O)=0.2, q(i,j) = 0.5, repeat the simulation N times (for N 
% at least 100) and plot the histograms of the states (open and closed) at 
% a few of the 1000 time steps (e.g. every 100). Report whether you 
% see convergence to the target distribution and approximately how quickly 
% it happens - this is called the burn-in time of the Markov chain simulation.

%% Part 2: Multiple ion channels
% For this section, always use the parameters: P(O) = 0.2 and q(i,j) = 0.5 
%
% Modify your code to increase the number of ion channels to a changeable 
% parameter M, with each channel opening and closing independently (this 
% only requires adding a loop). Save the total number of open and closed 
% channels at each time step.
%
% 2.1 Run the simulation at several values of M: 2, 10, 100 for 1000 time
% steps and report the mean number of open channels and the variance for
% each of the three cases. Is this what you expected? 
%
% 2.2 Repeat the multiple channel simulation N times (for N at least 100) and
% plot the histograms of the number of channels at several different time
% points (e.g. plot them at t=1, t=100, t=200, etc.
% Do they resemble any probability distribution you know?
%
% 2.3 We would like to find the burn-in time for converging to the
% stationary distribution. To do this, plot the
% *average* number of open/closed channels over the first 100 or so time
% steps. How quickly does this number reach an equilibrium value?

%
%% Part 3: 2 correlated ion channels
% Now let us investigate the situation when the probability of opening and 
% closing of a channel depends on the state of its neighbor. Let the 
% target joint probability distribution be:
%   P(O,O) = 0.6; P(C,C) = 0.2; P(O,C) = P(C,O) = 0.1
%
% Calculate the conditional probabilities for states of
% one channel given the state of its neighbor based on the joint distribution
% and use the Gibbs algorithm to compute a random sample. The Gibbs
% sampler does not accept or reject proposed states, but instead uses
% conditional probabilities to generate new states of the two random
% variables X and Y (the two ion channels): 
%
% Gibbs Sampler
%
% # Pick X0 and choose Y0 from the conditional distribution P (Y |X = X0) 
% # Choose X1 from P(X|Y = Y0) and Y1 from P(Y |X = X1)
% # Repeat N times
% 
% 3.1 Run for 1000 steps and plot the histograms for number of open channels 
% for the first 100 time steps and for the second 100 time steps. How
% different are they? Plot the histogram for 201 through 1000 time steps,
% and compare with the first two. Approximately how long does it take for
% the histrograms to converge?

%% Part 4: Multiple correlated ion channels
%
% Now we will use the Gibbs Sampler to simulate ion channel switching in
% the case that there are multiple correlated ion channels. The state of
% each ion channel depends on its right-hand and left-hand nearest neighbor: 
%
% P(O|2 open neighbors) = 0.4; P(O|1 open neighbor) = 0.2; P(O|0 open neighbors) = 0.1
% For the two ends of the ion channel array, the state of each ion channel
% can only depend on one nearest neighbor: 
% P (X = O|1 open neighbor) = 0.2 P (X = O|0 open neighbors) = 0.1
% The other probabilities are complementary.
%
% Gibbs Sampler with Multiple Correlated Channels
%
% # Generate an array of 100 ion channels in randomly-assigned states
% (0 or 1). 
% # Choose the state of the first ion channel in the array X(t,1) based on
% the state of its one neighbor X(t?1,2) from the previous time step. 
% # Choose the state of the second ion channel in the array X(t,2) based on 
% its neighbors X(t,1) and X(t?1,3), and do this for channels 2 through 99.
% # Select X(t,100) based on the state of its one neighbor X(t,99). 
% # Repeat N times.
%
% 4.1 Run for 10000 steps and plot the histogram for the number of channels 
% being open over all except for the first 100 time steps. 
% Approximately what number of open channels (total ion current) is
% expected? 
% 
% 4.2 Change your conditional probabilities to 
% P(O|2 open neighbors) = 0.2; P(O|1 open neighbor) = 0.5; 
% P(O|0 open neighbors) = 0.7 and re-run the simulation. Compare this 
% histogram to the one you obtained previously and estimate the expected
% current through the membrane.  

