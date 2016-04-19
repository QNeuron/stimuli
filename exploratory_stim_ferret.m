%% Reveil matin - exploratory stim
% Total duration 5 sec
% 2 sec noise burst
% 2 sec DRC
% 2 sec sweeps
% Music ?
% Voc ?
% Tones ?

sr = tdt100k;
finalSig = [];

%% White Noise
duration = 0.2; % second

totDur = 0;
sig = [];
while totDur < 2 %second
    y = wgn(ceil(duration*sr),1,1)'; % Generate white noise
    durSil = randi([50 300],1)*0.001; % second
    silence = ceil(durSil * sr);
    sig = [sig y./max(y) zeros(1,silence)];
    totDur = length(sig)/sr;
end

finalSig = [finalSig sig];

%% DRC (500Hz - 32kHz)

stim.f_s = tdt100k; % TDT 200K
stim.freq.min = 500;
stim.freq.spacing = 1/4; % octaves
stim.freq.step = 2^(stim.freq.spacing);
stim.freq.n = 6*(1/stim.freq.spacing)+1;

% get frequencies of tones
stim.freq.multipliers = 0:(stim.freq.n-1);
stim.freq.freqs = stim.freq.min*stim.freq.step.^stim.freq.multipliers;
stim.freqs = stim.freq.freqs;
fprintf('Frequency range %0.3f kHz to %0.3f kHz\n', min(stim.freqs/1000), max(stim.freqs/1000));

stim.ramp_duration = 15/1000; %seconds

stim.duration = 2; % sec
chord_durations = 100/1000; % msec

n_chords = stim.duration./chord_durations;

stim.meanlevel = 40;
% ranges = [10 20 30 40 99];
ranges = 40; % level-(range/2):level+(range/2)
switching_ranges = [20 40];

nToken = 25; % Number of random iterations

% make DRCs one at a time
stim.range = ranges;

stim.chord_duration = chord_durations;
stim.n_chords = n_chords;

for token = 1:nToken
    rand('seed', 110876+token*997);
    stim.token = token;
    fprintf('Range %d, chord duration %d, token %d', stim.range, stim.chord_duration*1000, stim.token);
    
    % make grid of levels
    stim.levels = (rand(stim.freq.n, stim.n_chords)-0.5)*stim.range+stim.meanlevel;
    
    % make waveform
    stim.drc = levels2drc(stim.f_s, stim.freqs, stim.levels, ...
        stim.chord_duration, stim.ramp_duration);
end

finalSig = [finalSig stim.drc.snd./max(stim.drc.snd)];

%% Sweep simple
freqs = [500 10000 5000 30000; ...
    5000 30000 500 10000];
nStim = size(freqs,2);

% length of a sweep in one direction
sweepDur = 0.4; % second
silDur = 0.2; % second

% number of samples of a sweep in one sweep direction
N = round(sweepDur * sr);
Nsil = round(silDur * sr);

% instantaneous frequency at each point in time ...
% ... first increasing for N points, then decreasing for N points
inst_f = zeros(nStim*N+nStim*Nsil,1);
c = 1;
for i = 1:nStim
    inst_f(c:c+N-1) = linspace(freqs(1,i), freqs(2,i), N);
    c = c+N+Nsil;
end
phi = 2 * pi * cumsum(inst_f) / sr;
sweep = sin(phi);

finalSig = [finalSig sweep'];
% pl = audioplayer(finalSig,sr);
% pl.play;


%% Sweeps continuous

% % parameters of linear up/down-sweep
% freq1 = 1000;
% freq2 = 30000;
% fs = sr;
% 
% % length of a sweep in one direction
% sweepDur = 3;
% 
% % number of samples of a sweep in one sweep direction
% N = round(sweepDur * fs);
% 
% % instantaneous frequency at each point in time ...
% % ... first increasing for N points, then decreasing for N points
% % inst_f = [linspace(freq1, freq2, N) linspace(freq2, freq1, N)];
% % inst_f = randi([500 10000],[1 ceil(N)]);
% % inst_f = [];
% % inst_f(1) = 500;
% % for i = 2:N,
% %     if rand>0.5
% %     inst_f(i) = inst_f(i-1) + randi([1 50],1);
% %     else
% %         inst_f(i) = inst_f(i-1) - randi([1 50],1);
% %     end
% % end
% 
% % Markov chain 2 - chain on freq
% Freqs = freq1:5:freq2;
% TRANS = diag(ones(length(Freqs),1)*0.5);
% pList = [0.2 0.1 0.1 0.05 0.05];
% for i = 1:length(pList),
%     TRANS = TRANS + diag(ones(length(Freqs)-i,1)*pList(i),i);
%     TRANS = TRANS + diag(ones(length(Freqs)-i,1)*pList(i),-i);
% end
% EMIS = ones(length(Freqs),1);
% [seq,states] = hmmgenerate(N-1,TRANS,EMIS);
% inst_f = Freqs(states);
% 
% 
% % since (in continuous time) instantaneous frequency is derivative of
% % phase, we have to compute the "integral" to get the phase for sin().
% phi = 2 * pi * cumsum(inst_f) / fs;
% sweep = sin(phi);
% % pl = audioplayer(sweep,sr);
% % pl.play;

