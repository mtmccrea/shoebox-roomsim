dbstop if error, clear all
%%% GENERATE TEST_SAMPLES_DRY
fs = 48e3;

%% SETUP

%data_path = ['..' filesep '..' filesep 'Resources' filesep 'data'];
hrtfs_path = '/Users/mccorml1/Documents/HRIRs_aalto2016/';
if ispc
    %samples_path = 'D:\Users\Archontis Politis\Google Drive\Audio\test_scenes'; % Windows Desktop
    %wav_path = 'C:\Users\Archontis Politis\Desktop\PHOCAHL_Tests';
elseif ismac
    %samples_path = '/Users/polarch/Google Drive/Audio/test_scenes'; % Macbook Air
    samples_path = '/Users/mccorml1/Documents/aalto_samples/test_scenes';
end
codepath = ['..' filesep 'Code'];
addpath(codepath)
% Write audio files ?
WRITE_FILES = 1;

% load HRTFs
sofa = loadSofaFile([hrtfs_path 'kemarhead_aalto2016.sofa']);
hrtf_dirs = sofa.SourcePosition(1:2,:).'; % degrees
hrtf_mtx = sofa.IR; 
itd = calcITDsfromHRIRs(hrtf_mtx,fs,0);  

% scenarios
sampleList = {
    '4SPEAKERS_maleEng_femEng_maleDan_femDan.wav',
    'BAND_shakers_bass_strings_drums.wav',
    'SOUNDSCAPE_voice_handclaps_fountain_piano.wav',
    [filesep 'ORCHESTRA' filesep 'bruckner_multichannel.wav']
};
% select which mono samples are used and at which order from the
% multichannel files
sourceIndex = {
    [1 2 4],    % maleEng, femAng, femDan
    [3 4 2 1],  % strings, drums, bass, shakers
    [3 1 2 4],  % fountain, voice, handclaps, piano
    [1:24]      % original order
};
nScenes = length(sampleList);

%%% SETUP GEOMETRY

% Put the sources in the first three scenarios at 1m (no attenuation)
% place sources at certain angles and distances from receiver
% (anticlockwise positive angles)
src_dirs = {
    [0 0; 70 0; -60 0],
    [60 10; 20 10; -20 10; -60 10],
    [110 0; 30 0; -45 0; -90 60]
};
src_dists = {
    ones(3,1),
    ones(4,1),
    ones(4,1)
};

% With the orchestra samples, let's use the absolute coords, since we have
% them
load([samples_path filesep 'ORCHESTRA' filesep 'orchestra_coords.mat'],'orchestraCenter2src_xyz')
nSrc = size(orchestraCenter2src_xyz,1);
% receiver with respect to orchestra center
orchestraCenter2rec_xyz = [-6 0 1.2];
% find distances and angles
rec2src_xyz = orchestraCenter2src_xyz - ones(nSrc,1)*orchestraCenter2rec_xyz;
src_dists{end+1} = sqrt(sum(rec2src_xyz.^2,2));
src_dirs{end+1} = unitCart2sph(rec2src_xyz)*180/pi;

%% RENDER SCENARIOS (BINAURAL)

for nc = 1:nScenes
    dists = src_dists{nc};
    dirs = src_dirs{nc};
    nSrc = size(dirs,1);    
    
    [src_sigs, fs_src] = audioread([samples_path filesep sampleList{nc}]);
    if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs
    src_sigs = src_sigs(:,sourceIndex{nc});
    lSig = size(src_sigs,1);
    
    %if all sources at same distance, ignore propagation delays and spherical attenuation    
    EQUIDISTANT = all(dists==dists(1));
    if EQUIDISTANT
        bin_sigs_dry{nc,1} = binauralize(src_sigs, dirs, hrtf_mtx, hrtf_dirs, fs, fs);
    else
        c = 343;        
        att = 1./dists;
        time = dists/c;
        time_samples = round(fs*time);
        src_rec_sigs = zeros(lSig+max(time_samples),nSrc);
        for ns=1:nSrc
            src_rec_sigs(time_samples(ns)+(1:lSig),ns) = att(ns)*src_sigs(:,ns);
        end
        bin_sigs_dry{nc,1} = binauralize(src_rec_sigs, dirs, hrtf_mtx, hrtf_dirs, fs, fs);
    end
    
    if WRITE_FILES
        bin_sigs_dry{nc,1} = 0.95* bin_sigs_dry{nc,1}/max(abs(bin_sigs_dry{nc,1}(:)));
        audiowrite([samples_path filesep sampleList{nc,1} '_Kemar217_bin_dry.wav'], bin_sigs_dry{nc,1}, fs, 'BitsPerSample', 24);
    end
end

%% ENCODE TO AMBISONICS AND RENDER
for nc = 1:nScenes
    dists = src_dists{nc};
    dirs = src_dirs{nc};
    nSrc = size(dirs,1);
    % find directions quantized to the HRTF grid, since that was used
    % for the binaural rendering, to avoid any directional bias
    [~, hrtf_closest_dirs_rad] = findClosestGridPoints(hrtf_dirs*pi/180, dirs*pi/180);
    hrtf_closest_dirs_deg = hrtf_closest_dirs_rad*180/pi;
    
    [src_sigs, fs_src] = audioread([samples_path filesep sampleList{nc}]);
    if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs
    src_sigs = src_sigs(:,sourceIndex{nc});
    lSig = size(src_sigs,1);
    
    %if all sources at same distance, ignore propagation delays and spherical attenuation    
    EQUIDISTANT = all(dists==dists(1));
    
    order = 5;
    if EQUIDISTANT
        Y = sqrt(4*pi)*getRSH(order, hrtf_closest_dirs_deg)';
        ambi_sigs_dry{nc} = src_sigs * Y;
    else
        c = 343;
        att = 1./dists;
        time = dists/c;
        time_samples = round(fs*time);
        src_rec_sigs = zeros(lSig+max(time_samples),nSrc);
        for ns=1:nSrc
            src_rec_sigs(time_samples(ns)+(1:lSig),ns) = att(ns)*src_sigs(:,ns);
        end
        Y = sqrt(4*pi)*getRSH(order, hrtf_closest_dirs_deg)';
        ambi_sigs_dry{nc} = src_rec_sigs * Y;
    end
     
    if WRITE_FILES
        ambi_sigs_dry{nc} = 0.95* ambi_sigs_dry{nc}/max(abs(ambi_sigs_dry{nc}(:)));
        audiowrite([samples_path filesep sampleList{nc} '_ambi_o' num2str(order) '_dry.wav'], ambi_sigs_dry{nc}, fs, 'BitsPerSample', 24);
    end
end

%%% GENERATE TEST_SAMPLES_REVERB

%% SETUP

data_path = ['..' filesep '..' filesep 'Resources' filesep 'data'];
hrtfs_path = '/Users/mccorml1/Documents/HRIRs_aalto2016/';
if ispc
    %samples_path = 'D:\Users\Archontis Politis\Google Drive\Audio\test_scenes'; % Windows Desktop
    %wav_path = 'C:\Users\Archontis Politis\Desktop\PHOCAHL_Tests';
elseif ismac
    samples_path = '/Users/mccorml1/Documents/aalto_samples/test_scenes';
    %samples_path = '/Users/polarch/Google Drive/Audio/test_scenes'; % Macbook Air
end
codepath = ['..' filesep 'Code'];
addpath(codepath)
% Write audio files ?
WRITE_FILES = 1;

% load HRTFs
sofa = loadSofaFile([hrtfs_path 'kemarhead_aalto2016.sofa']);
hrtf_dirs = sofa.SourcePosition(1:2,:).'; % degrees
hrtf_mtx = sofa.IR; 
itd = calcITDsfromHRIRs(hrtf_mtx,fs,0);  
grids{1} = hrtf_dirs*pi/180;
array_irs{1} = hrtf_mtx;

% scenarios
sampleList = {
    '4SPEAKERS_maleEng_femEng_maleDan_femDan.wav',
    'BAND_shakers_bass_strings_drums.wav',
    'SOUNDSCAPE_voice_handclaps_fountain_piano.wav',
    [filesep 'ORCHESTRA' filesep 'bruckner_multichannel.wav']
};
% select which mono samples are used and at which order from the
% multichannel files
sourceIndex = {
    [1 2 4],    % maleEng, femAng, femDan
    [3 4 2 1],  % strings, drums, bass, shakers
    [3 1 2 4],  % fountain, voice, handclaps, piano
    [1:24]      % 24ch orchestra - original order
};
nScenes = length(sampleList);

%%% SETUP GEOMETRY

% room dimensions
rooms = {
    [10.2  7.1  3.2], % office space
    [18.1 14.3 10.2], % small venue
    [31   23   19],   % large soundscape
    [35.7 19.8 17.4]  % some concert hall dims picked up for Musikvereinsaal
};

% desired RT60 (without taking air absorption into account), either per
% octave band or broadband
rt60 = {
    [0.80   0.64    0.56    0.48    0.40    0.32],  % office space
    [1.10   0.90    0.85    0.8     0.75    0.65],  % small venue
    0.7,                                            % large soundscape
    2.5                                             % concert hall
    %[2.7706    2.9330    2.9144    2.9332    2.6648    1.8987    1.1478] % Musikevereinsaal, measured by Lokki and co.
    };

% Put the sources in the first three scenarios at fixed distances,
% place sources at certain angles and distances from receiver
% (anticlockwise positive angles)
src_dirs = {
    [0 0; 70 0; -60 0],
    [60 10; 20 10; -20 10; -60 10],
    [110 0; 30 0; -45 0; -90 60]
};
src_dists = {
    ones(3,1),
    3*ones(4,1),
    6*ones(4,1)
};

% With the orchestra samples, let's use the absolute coords, since we have
% them
load([samples_path filesep 'ORCHESTRA' filesep 'orchestra_coords.mat'],'orchestraCenter2src_xyz')
nSrc = size(orchestraCenter2src_xyz,1);
% receiver with respect to orchestra center
orchestraCenter2rec_xyz = [-6 0 1.2];
% find distances and angles
rec2src_xyz = orchestraCenter2src_xyz - ones(nSrc,1)*orchestraCenter2rec_xyz;
src_dists{end+1} = sqrt(sum(rec2src_xyz.^2,2));
src_dirs{end+1} = unitCart2sph(rec2src_xyz)*180/pi;

% convert rec-to-source coords to room-based coords (measured from
% leftmost-bottom corner)
rec = {
    [3.0    3.4     1.7],
    [10.0   7.1     1.8],
    [13     11.5    9],
    [14     9.5     1.2]
    };
for nc = 1:nScenes
    % convert source directions from listener-centric to room-centric
    nSrc = size(src_dirs{nc},1);
    src_xyz = zeros(nSrc,3);
    [src_xyz(:,1), src_xyz(:,2), src_xyz(:,3)] = sph2cart(src_dirs{nc}(:,1)*pi/180, ...
        src_dirs{nc}(:,2)*pi/180, src_dists{nc});
    src_xyz(:,2) = -src_xyz(:,2);
    src{nc} = ones(nSrc,1)*rec{nc} + src_xyz;
    if any(any([src{nc}(:,1)>rooms{nc}(1), src{nc}(:,2)>rooms{nc}(2), src{nc}(:,3)>rooms{nc}(3)]))
        error('Source coordinates out of room boundaries')
    end
end

APPLY_AIR_ABSORPTION = true;

% SH orders for receivers
rec_orders = 5;

%% RENDER RIRs

for nc=1:nScenes
    nBands = length(rt60{nc});
    % lowest octave band
    band_centerfreqs = 125;
    for (nb=2:nBands) band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end
    % absorption for approximately achieving the RT60 above - row per band
    abs_wall = findAbsCoeffsFromRT(rooms{nc}, rt60{nc});
    % some room stats
    [~,  d_critical] = room_stats(rooms{nc}, abs_wall);
    nRec = size(rec{nc},1);
    nSrc = size(src_dirs{nc},1);
%     % plot scene
%     plot_scene(rooms{nc}, src{nc}, rec{nc});

    maxlim = max(rt60{nc}); % just stop if the echogram goes beyond that time ( or just set it to max(rt60) )
    for nb = 1:nBands
        if (rt60{nc}(nb)<maxlim) limits(nb) = rt60{nc}(nb);
        else limits(nb,1) = maxlim;
        end
    end
    
    %%% SH RENDERING
    % compute echograms for SH
    abs_echograms_sh = compute_echograms_sh(rooms{nc}, src{nc}, rec{nc}, abs_wall, limits, rec_orders);
    % render to RIRs
    temp_rirs = render_sh_rirs(abs_echograms_sh, band_centerfreqs, fs);
    if APPLY_AIR_ABSORPTION
        disp('Applying air absorption')
        nr=1;
        for ns=1:nSrc
            for nch = 1:size(temp_rirs,2)
                temp_rir = applyAirAbsorption(temp_rirs(:,nch,nr,ns), fs, 50);
                temp_rirs(:,nch,nr,ns) = temp_rir(1:size(temp_rirs,1));
            end
        end
    end
    sh_rirs{nc} = temp_rirs;
    
    
    %%% BINAURAL RENDERING
    % compute echograms for HRTF dirs
    abs_echograms_array = compute_echograms_arrays(rooms{nc}, src{nc}, rec{nc}, abs_wall, limits);
    % render to RIRs
    temp_rirs = render_array_rirs(abs_echograms_array, band_centerfreqs, fs, grids, array_irs);
    if APPLY_AIR_ABSORPTION
        disp('Applying air absorption') 
        for ns=1:nSrc
            for nch = 1:size(temp_rirs{1},2)
                temp_rir = applyAirAbsorption(temp_rirs{1}(:,nch,ns), fs, 50);
                temp_rir = temp_rir.';
                temp_rirs{1}(:,nch,ns) = temp_rir(1:size(temp_rirs{1},1));
            end
        end
    end
    bin_rirs{nc} = temp_rirs{1};
    
end

%% GENERATE SOUND SCENES
for nc=1:nScenes 
    
    [src_sigs, fs_src] = audioread([samples_path filesep sampleList{nc}]);
    if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs
    src_sigs = src_sigs(:,sourceIndex{nc});
    lSig = size(src_sigs,1);
    
    % get the source images and the final mix
    
    temp_bin_sigs_rev = apply_source_signals_arrays({bin_rirs{nc}}, src_sigs);
    bin_sigs_rev{nc} = temp_bin_sigs_rev{1};
    if WRITE_FILES
        bin_sigs_rev{nc} = 0.95* bin_sigs_rev{nc}/max(abs(bin_sigs_rev{nc}(:)));
        audiowrite([samples_path filesep sampleList{nc,1} '_Kemar217_bin_rev.wav'], bin_sigs_rev{nc}, fs, 'BitsPerSample', 24);
    end
    
    % get the source images and the final mix
    ambi_sigs_rev{nc} = apply_source_signals_sh(sh_rirs{nc}, src_sigs);

    if WRITE_FILES
        ambi_sigs_rev{nc} = 0.95* ambi_sigs_rev{nc}/max(abs(ambi_sigs_rev{nc}(:)));
        audiowrite([samples_path filesep sampleList{nc} '_ambi_o' num2str(order) '_rev.wav'], ambi_sigs_rev{nc}, fs, 'BitsPerSample', 24);
    end
end
