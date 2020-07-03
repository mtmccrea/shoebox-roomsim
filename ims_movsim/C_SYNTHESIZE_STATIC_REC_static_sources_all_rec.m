function [shrecsigs] = C_SYNTHESIZE_STATIC_REC_static_sources_all_rec(scene, src_sigs, sh_order, WRITE_TO_DISK)
%C_SYNTHESIZE_MOVING_REF Summary of this function goes here

if nargin<4
    WRITE_TO_DISK = 0;
end

hrtfs_path = '';

% load HRTFs
sofa = loadSofaFile([hrtfs_path 'kemarhead_aalto2016.sofa']);
hrtf_dirs{1} = sofa.SourcePosition(1:2,:).'; % degrees
hrtf_dirs{1} = hrtf_dirs{1} *pi/180;
hrtf_mtx{1} = sofa.IR; 
%itd = calcITDsfromHRIRs(hrtf_mtx{1},scene.fs,0);  

%% Generate room IRs at receiver points for interpolation

% limit the RIR by reflection order or by time-limit
nBands = length(scene.rt60);
maxlim = max(scene.rt60); % just cut if it's longer than that ( or set to max(rt60) )
for nb = 1:nBands
    if (scene.rt60(nb)<maxlim) limits(nb,1) = scene.rt60(nb);
    else limits(nb,1) = maxlim;
    end
end

% compute echograms
nCH = (sh_order+1)^2;
% change y coord for src/rec due to convention inside the IMS function
src2 = [scene.src(:,1) scene.room(2)-scene.src(:,2) scene.src(:,3)]; 
rec2 = [scene.rec(:,1) scene.room(2)-scene.rec(:,2) scene.rec(:,3)]; 

% SH echograms
abs_echograms = compute_echograms_sh(scene.room, src2, rec2, scene.abs_wall, limits, sh_order);
% array echograms from binaural reproduction
% abs_echograms_array = compute_echograms_arrays(scene.room, src2, rec2, scene.abs_wall, limits);

% render RIRs
sh_rirs = render_sh_rirs(abs_echograms, scene.abs_bands, scene.fs);
% bin_rirs = render_array_rirs(abs_echograms_array, scene.abs_bands, scene.fs, hrtf_dirs, hrtf_mtx);
nSrc = size(src2,1);
nMic = size(rec2,1);

if scene.AIR_ABS
    
    disp('Applying air absorption - SH')
%     sh_rirs_air = zeros(size(sh_rirs,1), nCH, nMic, nSrc);
    for ns=1:nSrc
        for nm=1:nMic
            for nc=1:nCH
                sh_rirs_air(:,nc,nm,ns) = applyAirAbsorption( ...
                    sh_rirs(:,nc,nm,ns), scene.fs);
            end
        end
    end
    shrecsigs = sh_rirs_air;
    
%     disp('Applying air absorption - Array')
%     bin_rirs_tmp = zeros(size(bin_rirs{1},1), size(bin_rirs{1},2), length(bin_rirs), size(bin_rirs{1},3));
%     for ns=1:nSrc
%         for nm=1:nMic
%             for nc = 1:size(bin_rirs{nm},2)
%                 temp_rir = applyAirAbsorption(bin_rirs{nm}(:,nc,ns), scene.fs);
%                 temp_rir = temp_rir.';
%                 bin_rirs_tmp(:,nc,nm,ns) = temp_rir(1:size(bin_rirs{nm},1));
%                %bin_rirs_tmp(:,nc,nm,ns) =  bin_rirs{nm}(:,nc,ns);
%             end
%         end
%     end 
%     bin_rirs = bin_rirs_tmp;
%     
%     clear bin_rirs_tmp;
end

% %% synthesize multiple receiver signals
% nRIR = size(scene.rec,1);


end

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

