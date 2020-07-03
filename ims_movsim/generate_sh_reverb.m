function [mixed_ir] = generate_sh_reverb(room_size,fc,rt60,src_pos,rec_pos,rec_order,fs)
%generate_sh_reverb Generate Spherical Harmonic domain impulse response.
% Combines Image Source Method for early part and exponential decay for
% late part of the impulse response.
%
% Syntax:  [mixed_ir] = generate_sh_reverb(room_size,fc,rt60,src_pos,rec_pos,rec_order,fs)
%
% Inputs:
%    room_size  : (3,)      - in m
%    fc         : (bands,)  - Center frequencies for reverb
%    rt60       : (bands,)  - RT_60 for each frequency band    
%    src_pos    : (nSRc,3)  - in m
%    rec_pos    : (3,)      - in m
%    rec_order  : int       - Spherical Harmonics order
%    fs         : int       - Sampling Frequency
%
% Outputs:
%    mixed_ir   : [t x nSH x nSrc]  - Impulse response
%
%
% Other m-files required: Spherical-Harmonic-Transform, shoebox-roomsim
%
% See also:
% ToDo: Check orthonormality limits, Check inverse factor
%
%   Chris Hold
%   Christoph.Hold@aalto.fi

PLOTFLAG = false;

%%
t = 0:1/fs:2*max(rt60)-1/fs;
nSmpl = length(t);
nSH = (rec_order+1)^2;
nBands = length(fc);


%% src rec conversions
[src_azi, src_ele, src_r] = cart2sph(src_pos(:, 1), ...
                                     src_pos(:, 2), ...
                                     src_pos(:, 3));
nSrc = size(src_pos,1);

% receiver position
nRec = size(rec_pos,1);

% % convert source directions from listener-centric to room-centric
[src_coords(:,1), src_coords(:,2), src_coords(:,3)] = sph2cart(src_azi, ...
    src_ele, src_r);
src = ones(nSrc,1)*rec_pos (1,:) + src_coords;
% check sources
for n=1:nSrc
    if (src(n,1)>room_size(1))||(src(n,2)>room_size(2))||(src(n,3)>room_size(3))||...
            (src(n,1)<0)||(src(n,2)<0)||(src(n,3)<0)
        error('Source coordinates out of room boundaries')
    end
end

%% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room_size, rt60);
if any(abs_wall>1)
    warning("Found absorption greater 1")
    abs_wall(abs_wall>1) = 1;
end

% critical distance for the room
[~,  d_critical] = room_stats(room_size, abs_wall);

l = room_size(1);
w = room_size(2);
h = room_size(3);

% total volume
V = l*w*h;

% total room wall area
Stot = 2*(l*w + l*h + w*h);

% Mixing time: Lindau, A., Kosanke, L., & Weinzierl, S. (2010). 
% Perceptual evaluation of physical predictors of the mixing time in binaural room impulse responses. 128th Audio Engineering Society Convention 2010, 3, 1418–1434.
t_mix = 20 * V/Stot + 12; % ms
disp(['Mixing time (ms): ' num2str(t_mix)])


%% rising ramp
nRamp = round(50/1000 * fs);  % 50ms
t_mix_a = ones(nSmpl, 1);
t_mix_smpl = round(t_mix/1000 *fs);
t_mix_a(1:t_mix_smpl-nRamp/2) = 0;
t_mix_a((t_mix_smpl-nRamp/2+1):(t_mix_smpl-nRamp/2+nRamp)) = linspace(0, 1, nRamp);

%% ISM Reverb
disp(num2str(nSrc) + " Sources, SH order: " + num2str(rec_order))

t_maxlim = (10*t_mix)/1000; % just stop if the echogram goes beyond that time ( or just set it to max(rt60) )
for nb = 1:nBands
    if (rt60(nb)<t_maxlim) limits(nb) = rt60(nb);
    else limits(nb,1) = t_maxlim;
    end
end
% ISM model has flipped y system.......
src_ism = ones(nSrc,1)*rec_pos (1,:) + ...
    [src_coords(:,1), -src_coords(:,2), src_coords(:,3)];

% compute echograms
[abs_echograms, ~,~] = compute_echograms_sh(room_size, src_ism, rec_pos, ...
                                            abs_wall, limits, rec_order);
% Render ISM SH IRs
ism_sh_rirs = render_sh_rirs(abs_echograms, fc, fs);
ism_sh_rirs = squeeze(ism_sh_rirs);
%assert(sum(imag(ism_sh_rirs), 'all') == 0)  % CFH: from absorption > 1

%% stochastic reverb
FLATTEN = 0;  % Produces ringin!
for iSrc = 1:nSrc
    stoch_rirs(:, :, iSrc) =  synthesizeNoiseReverb(nSH, fs, rt60, fc, FLATTEN);
end
assert(size(stoch_rirs, 1) <= nSmpl);
tmp = zeros(nSmpl, nSH, nSrc);
tmp(1:size(stoch_rirs, 1), :, :) = stoch_rirs;  % append zeros if necessary
stoch_rirs = tmp;
clear tmp

% randn is gaussian, normalize
stoch_rirs = stoch_rirs / max(max(abs(stoch_rirs(:, 1, :))));
% direct sound
[m1_ism, idx1_ism] = max(abs(squeeze(ism_sh_rirs(:, 1, :))));
% first reflection (offset by direct sound)
[m2_ism, idx2_ism] = max(abs(squeeze(ism_sh_rirs(idx1_ism+100:end, 1, :))));
idx_refl = idx1_ism + idx2_ism + 100;  % idx first reflection
m_refl = m2_ism;
assert(all(idx2_ism>100), "100 not enough")


%% Mix
mixed_ir = zeros(nSmpl,nSH,nSrc);

for iSrc = 1:nSrc
    % ISM part
    ism_part = ism_sh_rirs(:, :, iSrc) .* ...
        (1-t_mix_a(1:length(ism_sh_rirs)));  % apply slope/fade
    mixed_ir(1:length(ism_part),:,iSrc) = ism_part;
    % stochastic part
    stoch_part = cat(1, zeros(idx_refl(iSrc), nSH), stoch_rirs(:,:,iSrc));  % delay to first reflection
    stoch_part = stoch_part(1:nSmpl, :, :); % clip
    stoch_part = m_refl(iSrc) * stoch_part;  % match level to reflection
    stoch_part = (t_mix_a .* stoch_part);  % slope
    mixed_ir(:, :, iSrc) = mixed_ir(:, :, iSrc) + stoch_part;

end


%% PLTS
if PLOTFLAG
    figure(); hold on ; plot(ism_part(:, 1)); plot(stoch_part(:, 1));
    legend('ISM (p)', 'EXP (p)')
    title(iSrc)
    figure(); plot(t, mixed_ir(:, 1, iSrc))
end


%%
end  % //fnct


%% FROM HO-SIRR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rir_filt = synthesizeNoiseReverb(nCH, fs, t60, fc, FLATTEN)
%NOISEVERB Simulates a quick and dirty exponential decay reverb tail
%
% order:    HOA order
% fs:       sample rate
% t60:      reverberation times in different bands
% fc:       center frequencies of reverberation time bands (octave bands)
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

    if nargin<5, FLATTEN = 0; end
    
    % number of HOA channels
    nSH = nCH;
    % number of frequency bands
    nBands = length(t60);
    % decay constants
    alpha = 3*log(10)./t60;
    % length of RIR
    %lFilt = ceil(max(t60)*fs);
    t = (0:1/fs:2*max(t60)-1/fs)';  % go to 2*
    lFilt = length(t);
    % generate envelopes
    env = exp(-t*alpha);
    % generate RIRs
    rir = randn(lFilt, nSH, nBands);
    for k = 1:nBands
        rir(:, :, k) = rir(:,:,k).*(env(:,k)*ones(1,nSH));
    end
    % get filterbank IRs for each band
    filterOrder = 200;
    h_filt = filterbank(fc, filterOrder, fs);
    % filter rirs
    rir_filt = zeros(lFilt+ceil(filterOrder/2), nSH);
    for n = 1:nSH
        h_temp = [squeeze(rir(:,n,:)); zeros(ceil(filterOrder/2), nBands)];
        rir_filt(:, n) = sum(fftfilt(h_filt, h_temp), 2);
    end
                            
    if FLATTEN, rir_filt = equalizeMinphase(rir_filt); end
    
    rir_filt = rir_filt(filterOrder/2+1:end,:); % remove delay
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_filt = filterbank(fc, filterOrder, fs)
% fc:   the center frequencies of the bands
% Nord: order of hte FIR filter
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

    if length(fc) == 1
        h_filt = 1;

    elseif length(fc) == 2
        h_filt = zeros(filterOrder+1, 2);

        % lowpass
        f_ll = 2*fc(1)/sqrt(2);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fc(2)/sqrt(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, 2) = fir1(filterOrder, w_hh, 'high');

    else
        Nbands = length(fc);
        h_filt = zeros(filterOrder+1, Nbands);

        % lowpass
        f_ll = 2*fc(1)/sqrt(2);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fc(end)/sqrt(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, end) = fir1(filterOrder, w_hh, 'high');
        % bandpass
        for k = 2:Nbands-1
            fl = fc(k)/sqrt(2);
            fh = 2*fc(k)/sqrt(2);
            wl = fl/(fs/2);
            wh = fh/(fs/2);
            w = [wl wh];
            h_filt(:, k) = fir1(filterOrder, w, 'bandpass');
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rir_filt_flat = equalizeMinphase(rir_filt)
%MAKEFLATVERB Makes the decaying noise spectrally flat
%
%   Archontis Politis, 12/06/2018
%   archontis.politis@aalto.fi

Nrir = size(rir_filt,2);
for n=1:Nrir
    % equalise TDI by its minimum phase form to unity magnitude response
    tdi_f = fft(rir_filt(:,n));
    tdi_min_f = exp(conj(hilbert(log(abs(tdi_f)))));
    tdi_eq = real(ifft(tdi_f./tdi_min_f));
    rir_filt_flat(:,n) = tdi_eq;
end

end

