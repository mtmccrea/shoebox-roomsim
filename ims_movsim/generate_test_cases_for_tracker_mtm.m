clear all
close all

%%
path_base = '/Users/mccream1/Documents/Projects/6DOF/snd/';
path_src = [path_base 'src/'];
path_render = [path_base 'render/']; 

audiofiles = {'BAND_shakers_bass_strings_drums_short.wav'}; %, '4SPEAKERS_maleEng_femEng_maleDan_femDan.wav'};
audiofiles_shortnames = {'band'}; %, 'speech'};

%% SIMPLIFIED STATIC SOURCE/RECEIVER CASES
order = 2;
writeOutput = 0;
autoClose = 1;

% source positions (normalized to room dimensions)
clear src 

% src(2,:) = [0.8  0.7 0.35]; 
src(1,:) = [0.3  0.2 0.45]; 
% src(3,:) = [0.15 0.8 0.25]; 
% src(4,:) = [0.7  0.3 0.85]; 
Nsrc = size(src, 1);

% receiver positions (normalized to room dimensions)
rec(1,:) = [0.4 0.45 0.4];
rec(2,:) = [0.7 0.15 0.4];
Nrec = size(rec, 1);

for fi = 1:length(audiofiles)
    filepath = [ path_src  audiofiles{fi} ];  
    file_shortname = audiofiles_shortnames{fi};
    
    for si = Nsrc %1:4
        for roomi = {'anechoic'} %, 'medium', 'large'} 
            
            clear path_out
            path_out = [ path_render file_shortname '_' num2str(si) 'src_' ...
                num2str(Nrec) 'rec_static_o' num2str(order) '_' cell2mat(roomi) ];
            if ~exist(path_out, 'dir'), mkdir(path_out); end
            
            [scene, src_sigs] = A_SETUP_SCENE_static_src_rec( ...
                cell2mat(roomi), filepath, si, path_out, src, rec);
            
            % recsigs_sh [siglength,Nsh,Nrec,Nsrc]
            [recsigs_sh] = C_SYNTHESIZE_STATIC_REC_static_sources_all_rec(scene, src_sigs, order, 0);
            % sum sources across receivers
            recsigs_sh = sum(recsigs_sh, 4);
            
            [h_fig, h_mic, ~, ~] = plotScene(scene.room, scene.rec, scene.src, 1); % TOP VIEW ONLY for now
            hold on;
            set(gca, 'xlim', [-1 scene.room(1)+1], 'ylim', [-1 scene.room(2)+1])
            set(h_fig, 'Position', [65 292 907 647])
            
            if writeOutput
                % Requires R2020a or later
                exportgraphics(gcf, [scene.path_out filesep 'scene.pdf'], 'ContentType', 'image');
                
                % write each receiver signal to file
                for reci = 1:Nrec
                    audiowrite( ...
                        [scene.path_out filesep file_shortname '_o' num2str(order) '_ACN_N3D_' ...
                        cell2mat(roomi) '_' num2str(si) 'staticSrc_'  'staticRec' num2str(reci) '.wav' ], ...
                        recsigs_sh(:,:,reci), scene.fs, 'BitsPerSample', 24);
                end 

                % save metadata
                path_out = [path_out filesep 'metadata']; %#ok
                if ~exist(path_out, 'dir'), mkdir(path_out); end
                save([path_out filesep 'scene'], 'scene');
%                 save([path_out filesep 'trajectory'], 'trajectory');
            end
            if autoClose; pause(3), close(h_fig); end
        end

    end
end

%% Convert continuous signal to pressure-velocity components, set up windowing parameters

% convert to pressure velocity
M_sh2pv = beamWeightsPressureVelocity('real');
recsigs_pv = zeros(size(recsigs_sh,1), 4, Nrec);
for i = 1:Nrec
    recsigs_pv(:,:,i) = (M_sh2pv * recsigs_sh(:, 1:4, i)')';
end

% partition the signal into 100 sample buffers with 50% overlap
lBuf = 2048;                % buffer frame length
lHop = 1024;                % hop size
lOvlp = lBuf - lHop;        % length of overlap between frames
lSig = size(recsigs_pv,1);       % length of signal
Nbuf = ceil(lSig/lHop);     % number of frame buffers

pvSig_buffered = zeros(lBuf, 4, Nbuf, Nrec);

for reci=1:Nrec
    for chani=1:4
        pvSig_buffered(:,chani,:,reci) = buffer(recsigs_pv(:,chani,reci), lBuf, lOvlp);
    end
end

% Intensity Vector Calclulation of a continuous signals

% compute correlations of pressure-velocity for each partition, for
% short-time estimates of intensity vector
Ia_frm       = zeros(3, Nbuf, Nrec);
Ia_mag_frm   = zeros(1, Nbuf, Nrec);
Ia_err_frm   = zeros(1, Nbuf, Nrec);
diffuse_frm  = zeros(1, Nbuf, Nrec);
diffuse2_frm = zeros(1, Nbuf, Nrec);
DOA_cart_frm = zeros(3, Nbuf, Nrec);  % doa xyz
DOA_rad_frm  = zeros(2, Nbuf, Nrec);  % doa az, el
DOA_deg_frm  = zeros(2, Nbuf, Nrec);  % doa az, el
rms_frm      = zeros(1, Nbuf, Nrec);

eps = 1e-9;
% win = hann(lBuf);  % spottier results when windowing
win = ones(lBuf, 1); % rectangular window

for reci = 1:Nrec
    for bufi=1:Nbuf
        pvBlock = pvSig_buffered(:, :, bufi, reci);
        pvBlock = pvBlock .* win;
        p = pvBlock(:, 1);
        v = pvBlock(:, 2:4);

        % 3 Methods for active intensity vector
        % ~~~ Method 1: based on above instantaneous intensity vector calculation
        % Ia = sum(p .* v, 1);

        % ~~~ Method 2: similar to 1, but conjugate (to handle complex signal?)
        % see line 962 in TEST_SCRIPTS.m
        % NOTE: Intensity vector signals (actually they point to DoA instead of 
        % propagation as in the usual acoustic intensity convention)
        % Ia = sum(real((ones(3,1)*p) .* conj(v))'); 

        % ~~~ Method 3: correlations of SH signals
        % https://en.wikipedia.org/wiki/Covariance_matrix#Generalization_of_the_variance
        % the entries on the diagonal of the auto-covariance matrix 
        % K_XX are the variances of each element of the vector X
        % from getDiffuseness_IE.m
        %   sphCOV = (1/lSig) * (shsig*shsig');
        % correlations of PV signals
        % TODO: 0-mean signal: covariance = energy stat sig proc
        pvCOV = (1/sum(win)) * (pvBlock'*pvBlock); 
        Ia = real(pvCOV(2:4,1)); % [sum(p.*v_x), sum(p.*v_y), sum(p.*v_z)]

        % DOA and magnitude of intensity vector
        Ia_mag = norm(Ia);
        doa = Ia / Ia_mag;
        % normalized magnitude of intensity vector (experimental)
        Ia_mag_norm = Ia_mag / mean(p.^2); % experimental normalization (post doa calc)

        % energy and diffuseness estimate
        % Energy estimate TODO: why divide by 2?
        E = trace(pvCOV)/2 + eps;     % trace(pvCOV) = sum(sum(p.^2), sum(v_x.^2), sum(v_y.^2), sum(v_z.^2))
        diffuseness = 1 - (Ia_mag/E);
        tmp = dot(mean(v)', mean(v));
        tmp2 = mean(p.^2);
        diffuseness2 = 1 - (2 * Ia_mag / (tmp2 + tmp)); % eq. 4 in Leo's HO-SIRR 2020 paper, turns out this is the same as abs(Ia_mag_norm - 1)
        diffuseness2 = diffuseness2 * 0.5 + 0.5;         % hack to map [-1,1] > [0,1], TODO: why?

%         % Angular error (rad)
%         err = atan2(...
%             norm(cross(Ia, srcPos_fromReceiver)), ...
%             dot(Ia, srcPos_fromReceiver));

        % store results
        Ia_frm(:,bufi, reci)       = Ia;
        Ia_mag_frm(bufi, reci)     = Ia_mag_norm; % store normalized magnitude
%         Ia_err_frm(bufi, reci)     = rad2deg(err);
        diffuse_frm(bufi, reci)    = diffuseness;
        diffuse2_frm(bufi, reci)   = diffuseness2;
        DOA_cart_frm(:,bufi, reci) = doa;
        [DOA_rad_frm(1,bufi, reci), DOA_rad_frm(2,bufi)] = cart2sph(doa(1), doa(2), doa(3));
        rms_frm(bufi, reci)        = mag2db(sqrt(mean(p.^2)));
    end
end

DOA_deg_frm = rad2deg(DOA_rad_frm); % az, el

scene.doa = DOA_cart_frm;

% NOTE: See analyzeDecoder.m for proper use of sph2cart and angular error (angle
% between two vectors)

%% Insert DOAs into plot
[h_fig, h_mic, ~, ~] = plotScene(scene.room, scene.rec, scene.src, 1); % TOP VIEW ONLY for now
hold on;
set(gca, 'xlim', [-1 scene.room(1)+1], 'ylim', [-1 scene.room(2)+1])
set(h_fig, 'Position', [65 292 907 647])

for bufi = 1:Nbuf
    for i=1:Nrec
        % for top view
        stpt = scene.rec(i,1:2);
        endpt = stpt + scene.doa(1:2, bufi, i)'; % [[x,y,z],Nbuf,Nrec]
        plot([stpt(1) endpt(1)], [stpt(2) endpt(2)], '-');
    end
end
hold off;

%% Plot
% plot azim/elev with soundfile
% buf_stsamps = (0:nBuf-1) * (lHop+1) + 1;
buf_stsamps = (0:Nbuf-1) * (lHop) + 1;
buf_times = buf_stsamps / fs;

winheight = 0.1;

idx_odd  = 1:2:Nbuf;
idx_even = 2:2:Nbuf;

tmp1 = buf_stsamps(idx_odd); % starting points of windows
tmp2 = tmp1 + lBuf;          % ending points of windows
tmp1 = [tmp1; tmp2];
tmp1 = tmp1(:)';             % interleve: [start, end, start end, ...]
tmp1 = [tmp1; tmp1];         % stutter
winIdx_odd = tmp1(:)';       % and interleve again
winVal_odd = repmat([0 1 1 0], 1, ceil(numel(winIdx_odd)/4));
winVal_odd = winVal_odd(1:numel(winIdx_odd)) * winheight;

tmp1 = buf_stsamps(idx_even);
tmp2 = tmp1 + lBuf;
tmp1 = [tmp1; tmp2];
tmp1 = tmp1(:)';             % interleve: [start, end, start end, ...]
tmp1 = [tmp1; tmp1];         % stutter
winIdx_even = tmp1(:)';
winVal_even = repmat([0 1 1 0], 1, ceil(numel(winIdx_even)/4));
winVal_even = winVal_even(1:numel(winIdx_even)) * -winheight;

% mask to select above error threshold
doaerrthresh = 12; % error tolerance (deg)
rmsthresh = -40;
diffKeep_idx = Ia_err_frm < doaerrthresh;
rmsReject_idx = rms_frm < rmsthresh;

%TODO: circular mean
% CME from Tervo, Sakari. "Direction Estimation Based on Sound Intensity Vectors"
% angle(mean(exp(j*doas))))

figure;
maxlen = size(pvSig(:,1), 1) / fs;

% source signal
Nplt = 4;
p1 = subplot(Nplt,1,1); hold on;
plot((1:size(src_sig_mono(:,1), 1)) / fs, src_sig_mono/max(src_sig_mono) + 2, ...
    '-','Color', [1 1 1]*0.76, 'Marker', 'none', 'DisplayName', 'source'); % src mono signal
plot((1:size(src_sig_mono(:,1), 1)) / fs + 0.0163, src_sig_mono/max(src_sig_mono) + 2, ...
    '-','Color', hsv2rgb([0.8 0.8, 0.9]), 'Marker', 'none', 'DisplayName', 'source'); % src mono signal, delayed
stem((1:size(src_sig_sh(:,1), 1)) / fs, src_sig_sh(:,1)/max(src_sig_sh(:,1)), ...
   '-b', 'MarkerSize', 1.8, 'DisplayName', 'reverbed'); % omni from receiver signal
plot(winIdx_odd / fs,  winVal_odd,  'DisplayName', 'window', 'Selected','off', 'Visible', 'off'); % 'HandleVisibility', 'off', 
plot(winIdx_even / fs, winVal_even, 'DisplayName', 'window', 'Selected','off', 'Visible', 'off');
l = legend; hold off;
l.ItemHitFcn = @actionHighlight; % connect line visibility callback

% intensity vector AZIM, ELEV, ERROR
p2 = subplot(Nplt,1,2); hold on;
% title('Intensity AZIM, ELEV')
plot(buf_times, DOA_deg_frm(1,:), '.', 'Color', hsv2rgb([0 1 0.5]), 'DisplayName', 'az obser');
yline(src_dirs_deg(1), '--', 'Color', hsv2rgb([0.1 1 1]), 'DisplayName', sprintf("az enc %.1f", src_dirs_deg(1)));
plot(buf_times, DOA_deg_frm(2,:), '.', 'Color', hsv2rgb([0.4 1 1]), 'DisplayName', 'el obser'); hold on;
plot(buf_times, Ia_err_frm, 'r', 'DisplayName', 'I_a error (deg)');
yline(src_dirs_deg(2), '--', 'Color', hsv2rgb([0.5 1 0.5]), 'DisplayName',  sprintf("el enc %.1f", src_dirs_deg(2)));
legend; hold off;

% intensity vector MAG
p3 = subplot(Nplt,1,3); hold on;
plot((1:size(src_sig_sh(:,1), 1)) / fs, src_sig_sh/(2*max(src_sig_sh)) - 0.2, ...
    '-','Color', [1 1 1]*0.85, 'Marker', 'none', 'DisplayName', 'source del comp'); % src mono signal, delayed
yline(1, 'b');
yline(0, 'b');
yline(doaerrthresh / 180, 'r--');
% l1 = plot(buf_times, Ia_mag_frm, '-d', 'MarkerSize', 3, 'Color', '#7E2F8E', 'DisplayName', '||I_a|| / p^2');
% l2 = plot(buf_times, Ia_err_frm / 180, 'r', 'DisplayName', 'I_a error (norm)');
hsv = [0.4 0.3 0.8];
l1 = plot(buf_times, Ia_mag_frm, '-d', 'MarkerSize', 1.5, 'Color', hsv2rgb(hsv), 'DisplayName', '||I_a|| / p^2');

l2 = plot(buf_times, abs(diffuse2_frm-1), '-c', 'MarkerSize', 1.5, 'DisplayName', 'abs(\Psi_2-1)');

hsv(2) = 1; % full saturation for selected
plot(buf_times(diffKeep_idx), Ia_mag_frm(diffKeep_idx), 'd', ...
    'MarkerSize', 3, 'MarkerFaceColor',  hsv2rgb(hsv), 'MarkerEdgeColor', hsv2rgb(hsv));
plot(buf_times(rmsReject_idx), Ia_mag_frm(rmsReject_idx), 'rx', 'MarkerSize', 3);
l3 = plot(buf_times, Ia_err_frm / 180, 'r', 'DisplayName', 'I_a error (norm)');
legend([l1, l2, l3]); hold off;
% ylim([0 1] + ([-1 1] * 0.2))

% intensity vector DIFFUSENESS
p4 = subplot(Nplt,1,4); hold on;
% plot((1:size(src_sig_sh(:,1), 1)) / fs + 0.0163, src_sig_sh/(2*max(src_sig_sh)) - 0.2, ...
plot((1:size(src_sig_sh(:,1), 1)) / fs, src_sig_sh/(2*max(src_sig_sh)) - 0.2, ...
    '-','Color', [1 1 1]*0.85, 'Marker', 'none', 'DisplayName', 'source del comp'); % src mono signal, delayed
yline(1, 'b');
yline(0, 'b');
yline(doaerrthresh / 180, 'r--');
hsv = [0.6 0.5 0.8];
l1 = plot(buf_times, diffuse_frm, '-o', 'Color', hsv2rgb(hsv), 'MarkerSize', 1.5, 'DisplayName', '\Psi');
l2 = plot(buf_times, diffuse2_frm, 'Color', 'g', 'MarkerSize', 1.5, 'DisplayName', '\Psi-2');
hsv(2) = 1; % full saturation for selected
plot(buf_times(diffKeep_idx), diffuse_frm(diffKeep_idx), 'o', ...
    'MarkerSize', 3, 'MarkerFaceColor', hsv2rgb(hsv), 'MarkerEdgeColor', hsv2rgb(hsv));
plot(buf_times(rmsReject_idx), diffuse_frm(rmsReject_idx), 'rx', 'MarkerSize', 3);
l3 = plot(buf_times, Ia_err_frm / 180, 'r', 'DisplayName', 'I_a error (norm)');
legend([l1, l2, l3]); hold off;
% ylim([0 1] + ([-1 1] * 0.2))

% % intensity vector ERROR
% p4 = subplot(Nplt,1,4); hold on;
% % title('Intensity vector Error')
% plot(buf_times, Ia_err_frm, '-o', 'MarkerSize', 3, 'DisplayName', 'I_a error (deg)');
% legend; hold off;

linkaxes([p1 p2 p3 p4], 'x')
p1.XLim = [0 maxlen];

sgtitle(sprintf('Intensity metrics below %d˚ error threshold', doaerrthresh))
clear l1 l2
