function spikeTrain = binario(MUPulses,totalDur)

spikeTrain=zeros(length(MUPulses),totalDur);

for mu=1:length(MUPulses)
    
    spikeTrain(mu,MUPulses{mu})=1;

end

end

% A Cumulative Spike Train (CST) is the summation of individual motor unit spike trains identified through signal decomposition. 
% It serves as a representation of the collective activity of a motor unit pool and is often used to estimate the common synaptic input 
% (the "neural drive") to a muscle.
% First, you need to binarize the MUPulses file in a matrix so that you have the firings (see file binario attached above, it does exactly that).
% Then, after you have the firings you basically sum them:
% And you have a vector of the CST. Now from this CST, you have multiple options:
% 1) One option is to smooth it with a Hanning window:
% The smoothed CST is often highly correlated with the force or torque produced by the muscle, explaining up to 60-70% of the variance.
% Smooth the CST to estimate neural drive (e.g., 400 ms Hanning window):
% 2) another option is to do a window, sum all the firings in that window, and divide the sum by the number of active motor units 
% and divide again by the time window, so you have an average number of motor unit firing per second.

