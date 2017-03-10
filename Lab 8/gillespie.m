function gillespie(c, x, tmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function solves the Gillespie algorithm. It requires the rates and
% the initial states, along with a maximum time. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0; %sets initial time 
N= length(x);

while (t<tmax)
    r1 = rand;
    r2 = rand; 
    tau = -log(r1)/A;
end
end