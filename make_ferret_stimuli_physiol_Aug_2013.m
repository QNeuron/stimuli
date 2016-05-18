%% Stim generation
% this script generates stimuli for ferret behavioral pitch discrimination
% experiments
%
% Jan 7 2013
%
% Added a way to adjust RMS level for Benware. Quentin 2015.

close all
clear

sr = 97656;
stim_dur_ms = 200;

lowcutoff=200;
midcutoff=4000;
highcutoff=10000;

noise_l_co = 200;
noise_h_co = highcutoff;

stepsperoct=2;
numoct=2;
f0 = round(octavesteps(250,stepsperoct,stepsperoct*numoct+1));

rms_scale_factor = 100;

tone_adjustment=1;


% Amplitude adjustment params

wantedDB = 80;
BenwaredBrms1 = 94;
P0 = 1;
fileFormat = 'f32';

% Stage 2 - tones with high harmonics
%
% 3000-7200 Hz
% tone 1: 300 Hz F0, harmonics 10-24; tone 2: 600 Hz, harmonics 5-12

l_co = midcutoff;
u_co = highcutoff;

tone_hh = [];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        if size(tone_hh,1)<ii
            tone_hh(ii,:) = harm;
        else
            tone_hh(ii,:) = zadd(tone_hh(ii,:),harm);
        end
    end
    if ii==1,
        tone_hh(ii,:) = tone_hh(ii,:)/rms_scale_factor;
    else
        tone_hh(ii,:) = tone_hh(ii,:)/rms(tone_hh(ii,:))*rms(tone_hh(1,:));
    end
end

dp_noise = pnoise(stim_dur_ms,noise_l_co,noise_h_co,-30,0,sr);

% ERBs for ferret: 180 Hz at 300 Hz, 225 Hz at 500 Hz, 250 Hz at 600 Hz, 300 Hz at 1000 Hz
%cf = 300;
%erb_width = 180;
%cf=600; erb_width=250;
%calibrating with level at a higher F0 is more conservative as filter there is
%apparently narrower ==> less noise power
cf = 1000;
erb_width = 300;
erb_low = (-erb_width + sqrt(erb_width^2 + 4*cf^2))/2;
erb_high = erb_low+erb_width;

dp_noise_fft = fft(dp_noise);

%signal is odd length
nfreqs = (length(dp_noise)-1)/2;
max_freq = sr*(length(dp_noise)-1)/2/length(dp_noise); %max freq is just under nyquist
freqs = [0:max_freq/nfreqs:max_freq];
neg_freqs = fliplr(freqs(2:end));

[temp, low_bin] = min(abs(freqs-erb_low));
[temp, high_bin] = min(abs(freqs-erb_high));

[temp, low_bin_neg] = min(abs(neg_freqs-erb_low));
[temp, high_bin_neg] = min(abs(neg_freqs-erb_high));
low_bin_neg = low_bin_neg+length(freqs);
high_bin_neg = high_bin_neg+length(freqs);

dp_noise_fft_filt = dp_noise_fft;
dp_noise_fft_filt([1:low_bin-1 high_bin+1:high_bin_neg-1 low_bin_neg+1:length(dp_noise)]) = 0;
dp_noise_filt = ifft(dp_noise_fft_filt);

%determine factor by which noise band must be multiplied to be 5 dB below
%primaries, i.e. 10 dB above hypothetical DPs which are supposed to be 15
%dB below primaries (based on Pressnitzer and Patterson data)
%
% P&P data show that DPs grow with number of primaries
% DP amplitude increases ~3dB for each doubling of the number of primaries
% in their experiment primaries were 60 db SPL and with 17 primaries (the
% most they tried), the DP was ~45 db SPL.
% Our high-harmonic tones have 15 or 8 harmonics (low/high f0), so
% accounting for a DP that is 15 dB below the primaries seems conservative.

sample_harmonic = tone(f0(1)*10,stim_dur_ms,0,sr);
rms_primary = rms(sample_harmonic)/rms_scale_factor;
noise_adjustment = 1/rms(dp_noise_filt) * rms_primary * 10^(-5/20);
dp_noise_adjusted = dp_noise*noise_adjustment*tone_adjustment;

tone_high_harm_stim=[];
for ii=1:length(f0),
    tone_high_harm_stim(ii,:) = zadd(tone_hh(ii,:), dp_noise_adjusted);
end


% Stage 1 - full bandwidth tones
%
% 600-7200 Hz
% tone 1: 300 Hz F0, harmonics 2-24; tone 2: 600 Hz, harmonics 1-12

l_co = lowcutoff;
u_co = highcutoff;

tone_ah = [];
tone_all_harm_stim=[];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        if size(tone_ah,1)<ii
            tone_ah(ii,:) = harm;
        else
            tone_ah(ii,:) = zadd(tone_ah(ii,:),harm);
        end
    end
    tone_ah(ii,:) = tone_ah(ii,:)/rms(tone_ah(ii,:))*rms(tone_hh(1,:));
    tone_all_harm_stim(ii,:) = zadd(tone_ah(ii,:), dp_noise_adjusted);
end



% Stage 3 - tones with lower harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5

l_co = lowcutoff;
u_co = midcutoff;

tone_lh = [];
tone_low_harm_stim=[];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        if size(tone_lh,1)<ii
            tone_lh(ii,:) = harm;
        else
            tone_lh(ii,:) = zadd(tone_lh(ii,:),harm);
        end
    end
    tone_lh(ii,:) = tone_lh(ii,:)/rms(tone_lh(ii,:))*rms(tone_hh(1,:));
    tone_low_harm_stim(ii,:) = zadd(tone_lh(ii,:), dp_noise_adjusted);
end


% Stage 4 - ALT phase tones with higher harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5
l_co = midcutoff;
u_co = highcutoff;

tone_hh_alt = [];
tone_high_harm_alt_stim=[];
for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        if mod(h,2)==1,
            harm = tone(f0(ii)*h,stim_dur_ms,0,sr);
        else
            harm = tone(f0(ii)*h,stim_dur_ms,pi,sr);
        end
        if size(tone_hh_alt,1)<ii
            tone_hh_alt(ii,:) = harm;
        else
            tone_hh_alt(ii,:) = zadd(tone_hh_alt(ii,:),harm);
        end
    end
    tone_hh_alt(ii,:) = tone_hh_alt(ii,:)/rms(tone_hh_alt(ii,:))*rms(tone_hh(1,:));
    tone_high_harm_alt_stim(ii,:) = zadd(tone_hh_alt(ii,:), dp_noise_adjusted);
end


% Stage 5 - RAND phase tones with higher harmonics
%
% 600-3000 Hz
% tone 1: 300 Hz F0, harmonics 2-10; tone 2: 600 Hz, harmonics 1-5
l_co = midcutoff;
u_co = highcutoff;

tone_hh_rand = [];
tone_high_harm_rand_stim=[];

for ii=1:length(f0),
    for h = l_co/f0(ii) : u_co/f0(ii)
        thisphase=rand*2*pi;
    %     thisphase=rand;
        harm = tone(f0(ii)*h,stim_dur_ms,thisphase,sr);
        if size(tone_hh_rand,1)<ii
            tone_hh_rand(ii,:) = harm;
        else
            tone_hh_rand(ii,:) = zadd(tone_hh_rand(ii,:),harm);
        end
    end
    tone_hh_rand(ii,:) = tone_hh_rand(ii,:)/rms(tone_hh_rand(ii,:))*rms(tone_hh(1,:));
    tone_high_harm_rand_stim(ii,:) = zadd(tone_hh_rand(ii,:), dp_noise_adjusted);
end

% Just pure tones
tone_pure = [];
tone_pure_stim=[];
for ii=1:length(f0),
    tone_pure(ii,:) = tone(f0(ii),stim_dur_ms,0,sr);
    tone_pure(ii,:) = tone_pure(ii,:)/rms(tone_pure(ii,:))*rms(tone_hh(1,:));
    tone_pure_stim(ii,:) = zadd(tone_pure(ii,:), dp_noise_adjusted);
end

gatelength=5; %ramp time in msec
for ii=1:length(f0)
    tone_high_harm_rand_stim(ii,:)=linear_envelope(tone_high_harm_rand_stim(ii,:),gatelength,sr);
    tone_high_harm_alt_stim(ii,:)=linear_envelope(tone_high_harm_alt_stim(ii,:),gatelength,sr);
    tone_low_harm_stim(ii,:)=linear_envelope(tone_low_harm_stim(ii,:),gatelength,sr);
    tone_high_harm_stim(ii,:)=linear_envelope(tone_high_harm_stim(ii,:),gatelength,sr);
    tone_all_harm_stim(ii,:)=linear_envelope(tone_all_harm_stim(ii,:),gatelength,sr);
    tone_pure_stim(ii,:)=linear_envelope(tone_pure_stim(ii,:),gatelength,sr);
    tone_pure(ii,:)=linear_envelope(tone_pure(ii,:),gatelength,sr);
    tone_ah(ii,:)=linear_envelope(tone_ah(ii,:),gatelength,sr);
end

% Adjust the amplitude to match the wantedDB var;

for i = 1:length(f0),
        
    cRMS = rms(tone_pure_stim(i,:));
%     currentDB = BenwaredBrms1 + 20 * log10(cRMS);
    wRMS = 10^((wantedDB - BenwaredBrms1)/20);
    tone_pure_stim(i,:) = tone_pure_stim(i,:) .* (wRMS/cRMS);
    
    cRMS = rms(tone_all_harm_stim(i,:));
    wRMS = 10^((wantedDB - BenwaredBrms1)/20);
    tone_all_harm_stim(i,:) = tone_all_harm_stim(i,:) .* (wRMS/cRMS);
    
    cRMS = rms(tone_high_harm_stim(i,:));
    wRMS = 10^((wantedDB - BenwaredBrms1)/20);
    tone_high_harm_stim(i,:) = tone_high_harm_stim(i,:) .* (wRMS/cRMS);
    
    cRMS = rms(tone_low_harm_stim(i,:));
    wRMS = 10^((wantedDB - BenwaredBrms1)/20);
    tone_low_harm_stim(i,:) = tone_low_harm_stim(i,:) .* (wRMS/cRMS);
    
    cRMS = rms(tone_high_harm_alt_stim(i,:));
    wRMS = 10^((wantedDB - BenwaredBrms1)/20);
    tone_high_harm_alt_stim(i,:) = tone_high_harm_alt_stim(i,:) .* (wRMS/cRMS);
    
    cRMS = rms(tone_high_harm_rand_stim(i,:));
    wRMS = 10^((wantedDB - BenwaredBrms1)/20);
    tone_high_harm_rand_stim(i,:) = tone_high_harm_rand_stim(i,:) .* (wRMS/cRMS);
    
end

msg = 'Amplitude above 1 or bellow -1. WAV files will be clipped. Use fileFormat = f32 instead.';
if strcmpi(fileFormat,'wav');
    
    if max(tone_pure_stim(:))>1 || min(tone_pure_stim(:))<-1,
        warning();
        break;
    end
    if max(tone_all_harm_stim(:))>1 || min(tone_all_harm_stim(:))<-1,
        warning(msg);
        break;
    end
    if max(tone_high_harm_stim(:))>1 || min(tone_high_harm_stim(:))<-1,
        warning(msg);
        break;
    end
    if max(tone_low_harm_stim(:))>1 || min(tone_low_harm_stim(:))<-1,
        warning(msg);
        break;
    end
    if max(tone_high_harm_alt_stim(:))>1 || min(tone_high_harm_alt_stim(:))<-1,
        warning(msg);
        break;
    end
    if max(tone_high_harm_rand_stim(:))>1 || min(tone_high_harm_rand_stim(:))<-1,
        warning(msg);
        break;
    end
    
end
disp('done')

%% Plot Spectrum

% plotting stimuli
figpos=[1283 -233 1662 420];
xx=highcutoff/1000+1;
midstim=ceil(size(tone_all_harm_stim,1)/2);
xtype='linear';

figure %spectra
set(gcf,'position',figpos)
subplot(3,5,1);
pwelch(tone_all_harm_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title(f0(1))
subplot(3,5,2);
pwelch(tone_high_harm_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title('High harmonics')
subplot(3,5,3);
pwelch(tone_low_harm_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title('Low harmonics')
subplot(3,5,4);
pwelch(tone_high_harm_alt_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title('ALT phase')
subplot(3,5,5);
pwelch(tone_high_harm_rand_stim(1,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title('Random phase')

subplot(3,5,6);
pwelch(tone_all_harm_stim(midstim,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title(f0(midstim))
subplot(3,5,7);
pwelch(tone_high_harm_stim(midstim,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
subplot(3,5,8);
pwelch(tone_low_harm_stim(midstim,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
subplot(3,5,9);
pwelch(tone_high_harm_alt_stim(midstim,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
subplot(3,5,10);
pwelch(tone_high_harm_rand_stim(midstim,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);

subplot(3,5,11);
pwelch(tone_all_harm_stim(end,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
title(f0(end))
subplot(3,5,12);
pwelch(tone_high_harm_stim(end,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
subplot(3,5,13);
pwelch(tone_low_harm_stim(end,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
subplot(3,5,14);
pwelch(tone_high_harm_alt_stim(end,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);
subplot(3,5,15);
pwelch(tone_high_harm_rand_stim(end,:),[],[],[],sr);
set(gca,'YLim',[-100 -50],'XLim',[0 xx],'XScale',xtype);



% figure %spectra on a log scale
% set(gcf,'position',figpos)
% subplot(2,5,1);
% pwelch(tone1_all_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,2);
% pwelch(tone1_high_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,3);
% pwelch(tone1_low_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,4);
% pwelch(tone1_high_harm_alt_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,5);
% pwelch(tone1_high_harm_rand_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,6);
% pwelch(tone2_all_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,7);
% pwelch(tone2_high_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,8);
% pwelch(tone2_low_harm_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,9);
% pwelch(tone2_high_harm_alt_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
% subplot(2,5,10);
% pwelch(tone2_high_harm_rand_stim,[],[],[],sr);
% set(gca,'YLim',[-100 -50],'XLim',[.1 8],'XScale','log');
%% Plot waveform


figure %waveform
set(gcf,'position',figpos)
% plotdur=length(tone1_all_harm_stim);%duration to plot
plotdur=floor(.1*sr); %plot just first 'x' ms
subplot(2,5,1);
plot([1:plotdur]./sr,tone_all_harm_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2]);
title('All harmonics')
subplot(2,5,2);
plot([1:plotdur]./sr,tone_high_harm_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2]);
title('High harmonics')
subplot(2,5,3);
plot([1:plotdur]./sr,tone_low_harm_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2]);
title('Low harmonics')
subplot(2,5,4);
plot([1:plotdur]./sr,tone_high_harm_alt_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2]);
title('ALT phase')
subplot(2,5,5);
plot([1:plotdur]./sr,tone_high_harm_rand_stim(1,1:plotdur));
set(gca,'YLim',[-.2 .2]);
title('Random phase')
subplot(2,5,6);
plot([1:plotdur]./sr,tone_all_harm_stim(end,1:plotdur));
set(gca,'YLim',[-.2 .2]);
subplot(2,5,7);
plot([1:plotdur]./sr,tone_high_harm_stim(end,1:plotdur));
set(gca,'YLim',[-.2 .2]);
subplot(2,5,8);
plot([1:plotdur]./sr,tone_low_harm_stim(end,1:plotdur));
set(gca,'YLim',[-.2 .2]);
subplot(2,5,9);
plot([1:plotdur]./sr,tone_high_harm_alt_stim(end,1:plotdur));
set(gca,'YLim',[-.2 .2]);
subplot(2,5,10);
plot([1:plotdur]./sr,tone_high_harm_rand_stim(end,1:plotdur));
set(gca,'YLim',[-.2 .2]);

xlabel('Time (s)')
axis 'tight'


%% listening to phase manipulations

% fullstimHIGH=[tone_all_harm_stim(midstim,:) zeros(1,ceil(sr*.5)) tone_high_harm_stim(midstim,:) zeros(1,ceil(sr*.5))];
% fullstimLOW=[tone_all_harm_stim(midstim,:) zeros(1,ceil(sr*.5)) tone_low_harm_stim(midstim,:) zeros(1,ceil(sr*.5))];
% fullstimALT=[tone_high_harm_stim(midstim,:) zeros(1,ceil(sr*.5)) tone_high_harm_alt_stim(midstim,:) zeros(1,ceil(sr*.5))];
% fullstimRAND=[tone_high_harm_stim(midstim,:) zeros(1,ceil(sr*.5)) tone_high_harm_rand_stim(midstim,:) zeros(1,ceil(sr*.5))];
% 
% soundsc(repmat(fullstimHIGH,1,5),sr)
% 
% soundsc(repmat(fullstimLOW,1,5),sr)
% 
% soundsc(repmat(fullstimALT,1,5),sr)
% 
% soundsc(repmat(fullstimRAND,1,5),sr)

%% Save sounds
% close all
% save KerryPitchSoundsPhysiology2013.mat tone_pure_stim tone_all_harm_stim tone_high_harm_stim tone_low_harm_stim tone_high_harm_alt_stim tone_high_harm_rand_stim...
%     sr stim_dur_ms f0...
%     noise_l_co noise_h_co gatelength...
%     lowcutoff midcutoff highcutoff;

%% Save sound sequences for Benware
soundtypes=1:6; % 1 tone, 2 all harm, 3 high, 4 low, 5 alt, 6 rand 
reps=20;
alltrialf0s=[];
alltrialtypes=[];

for kk=1:reps
    for ii=1:length(f0)
        for jj=1:length(soundtypes)
            alltrialf0s=[alltrialf0s f0(ii)];
        end
        alltrialtypes=[alltrialtypes soundtypes];
    end
end

r=randperm(length(alltrialtypes));
alltrialf0s=alltrialf0s(r);
alltrialtypes=alltrialtypes(r);
acquisitiontime=length(alltrialtypes)/60 % acquisition total time, in minutes.
numfiles=length(alltrialtypes)/40; % Number of 40-sec files we need to make !!This has to be an integer

for ii=1:numfiles %for each file
    stimulus=[];
    stimInd=[1:40]+(ii-1)*40;
    trialf0s=alltrialf0s(stimInd);
    trialtypes=alltrialtypes(stimInd);
    for jj=1:length(trialf0s)
        switch trialtypes(jj)
            case 1
                stimulus=[stimulus tone_pure_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
            case 2
                stimulus=[stimulus tone_all_harm_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
            case 3
                stimulus=[stimulus tone_high_harm_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
            case 4
                stimulus=[stimulus tone_low_harm_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
            case 5
                stimulus=[stimulus tone_high_harm_alt_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
            case 6
                stimulus=[stimulus tone_high_harm_rand_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
        end
    end
    filename=['KerryPitchSounds2013_' num2str(ii)];
    switch fileFormat
        case 'wav'
            wavwrite(stimulus,sr,[filename '.wav']);
        case 'f32'
            fid = fopen([filename '.f32'], 'w');
            fwrite(fid,stimulus, 'float32');
            fclose(fid);
        otherwise
            error('Invalid file format.');
    end
    save([filename '.mat'],'stimulus','trialf0s','trialtypes','soundtypes','reps','alltrialf0s','alltrialtypes');
end
%% Sanity check
% sanity check 1
figure(1); clf
plot((1:length(stimulus))./sr,stimulus);
figure(2); clf
plot(trialf0s,'ok');

% sanity check 2
for ii=1
    stimulus=[];
    stimInd=[1:40]+(ii-1)*40;
    trialf0s=alltrialf0s(stimInd);
    trialtypes=alltrialtypes(stimInd);
    for jj=1:length(trialf0s)
        switch trialtypes(jj)
            case 1
                stimulus=[tone_pure_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
                plot((1:length(stimulus))./sr,stimulus)
                title([num2str(trialf0s(jj)) 'Hz, ' num2str(trialtypes(jj))])
                soundsc(stimulus,sr)
                pause()
            case 2
                stimulus=[tone_all_harm_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
                plot((1:length(stimulus))./sr,stimulus)
                title([num2str(trialf0s(jj)) 'Hz, ' num2str(trialtypes(jj))])
                soundsc(stimulus,sr)
                pause()
            case 3
                stimulus=[tone_high_harm_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
               
                plot((1:length(stimulus))./sr,stimulus)
                title([num2str(trialf0s(jj)) 'Hz, ' num2str(trialtypes(jj))])
                soundsc(stimulus,sr)
                pause()
            case 4
                stimulus=[tone_low_harm_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
                
                plot((1:length(stimulus))./sr,stimulus)
                title([num2str(trialf0s(jj)) 'Hz, ' num2str(trialtypes(jj))])
                soundsc(stimulus,sr)
                pause()
            case 5
                stimulus=[tone_high_harm_alt_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
                
                plot((1:length(stimulus))./sr,stimulus)
                title([num2str(trialf0s(jj)) 'Hz, ' num2str(trialtypes(jj))])
                soundsc(stimulus,sr)
                pause()
            case 6
                stimulus=[tone_high_harm_rand_stim(find(f0==trialf0s(jj)),:) zeros(1,floor((1000-stim_dur_ms)/1000*sr))];
                
                plot((1:length(stimulus))./sr,stimulus)
                title([num2str(trialf0s(jj)) 'Hz, ' num2str(trialtypes(jj))])
                soundsc(stimulus,sr)
                pause()
        end
    end
end

%% Stim spectral enveloppe
resolution = 500; % no oscillations above a freq > 1/resolution
stim = tone_high_harm_rand_stim;
% tone_pure_stim
% tone_all_harm_stim
% tone_high_harm_stim
% tone_low_harm_stim
% tone_high_harm_alt_stim
% tone_high_harm_rand_stim

for i = 1:size(stim,1),

    x = stim(i,:);
    
    disp(['F0 : ' num2str(f0(i)) ' RMS=' num2str(rms(x))])
    
    figure;
    title(['FO : ' num2str(f0(i))]);
    [y, f] = plot_fft(x,sr,20000);
    [yymax,yymin,maxi,argmaxi,mini,argmini]=enveloppe(y,resolution);
    hold on
    plot(f,yymax,'k');
    hold off
end

% Resolved harmonic (nerve model)

% All Tone
Plot_physiol_stim_resolvability(f0,lowcutoff,highcutoff)

% Low cut
Plot_physiol_stim_resolvability(f0,lowcutoff,midcutoff)

% High cut
Plot_physiol_stim_resolvability(f0,midcutoff, highcutoff)
