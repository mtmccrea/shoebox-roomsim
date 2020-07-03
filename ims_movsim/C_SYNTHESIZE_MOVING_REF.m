function [shmovsig, binmovsig] = C_SYNTHESIZE_MOVING_REF(scene, trajectory, src_sigs, sh_order, WRITE_TO_DISK)
%C_SYNTHESIZE_MOVING_REF Summary of this function goes here

if nargin<4
    WRITE_TO_DISK = 0;
end

%hrtfs_path = '/Users/mccorml1/Documents/HRIRs_aalto2016/';
hrtfs_path = '';

% load HRTFs
sofa = loadSofaFile([hrtfs_path 'kemarhead_aalto2016.sofa']);
hrtf_dirs{1} = sofa.SourcePosition(1:2,:).'; % degrees
hrtf_dirs{1} = hrtf_dirs{1} *pi/180;
hrtf_mtx{1} = sofa.IR; 
%itd = calcITDsfromHRIRs(hrtf_mtx{1},scene.fs,0);  
 

%% Generate room IRs at segment points for interpolation
% number of rir points
nRIR = size(trajectory.rir_pts,1);

% limit the RIR by reflection order or by time-limit
nBands = length(scene.rt60);
maxlim = max(scene.rt60); % just cut if it's longer than that ( or set to max(rt60) )
for nb = 1:nBands
    if (scene.rt60(nb)<maxlim) limits(nb) = scene.rt60(nb);
    else limits(nb,1) = maxlim;
    end
end
% compute echograms
nCH = (sh_order+1)^2;
src2 = [scene.src(:,1) scene.room(2)-scene.src(:,2) scene.src(:,3)]; % change y coord for src/rec due to convention inside the IMS function
rec2 = [trajectory.rir_pts(:,1) scene.room(2)-trajectory.rir_pts(:,2) trajectory.rir_pts(:,3)]; % change y coord for src/rec due to convention inside the IMS function
abs_echograms = compute_echograms_sh2(scene.room, src2, rec2, scene.abs_wall, limits, sh_order);
abs_echograms_array = compute_echograms_arrays(scene.room, src2, rec2, scene.abs_wall, limits);

% render RIRs
sh_rirs = render_sh_rirs(abs_echograms, scene.abs_bands, scene.fs);
bin_rirs = render_array_rirs_same_rec(abs_echograms_array, scene.abs_bands, scene.fs, hrtf_dirs, hrtf_mtx);


clear abs_echograms abs_echograms_array;
nSrc = size(src2,1);
nMic = size(rec2,1);
if scene.AIR_ABS
    disp('Applying air absorption - SH')
    for ns=1:nSrc
        for nm=1:nMic
            for nc=1:nCH
                sh_rirs_air(:,nc,nm,ns) = applyAirAbsorption(sh_rirs(:,nc,nm,ns),scene.fs);
            end
        end
    end
    sh_rirs = sh_rirs_air;
    clear sh_rirs_air;
    
    bin_rirs_tmp = zeros(size(bin_rirs{1},1), size(bin_rirs{1},2), length(bin_rirs), size(bin_rirs{1},3));
    disp('Applying air absorption - Array') 
    for ns=1:nSrc
        for nm=1:nMic
            for nc = 1:size(bin_rirs{nm},2)
                temp_rir = applyAirAbsorption(bin_rirs{nm}(:,nc,ns), scene.fs);
                temp_rir = temp_rir.';
                bin_rirs_tmp(:,nc,nm,ns) = temp_rir(1:size(bin_rirs{nm},1));
                
               %bin_rirs_tmp(:,nc,nm,ns) =  bin_rirs{nm}(:,nc,ns);
            end
        end
    end 
    
    bin_rirs = bin_rirs_tmp;
    clear bin_rirs_tmp;
end

%% synthesize moving receiver signal

addpath('stft_conv')

%duration of moving source signal
lMovSig = round(trajectory.time*scene.fs);
% replicate the source signals to fit in the moving one
src_mov_sigs = zeros(lMovSig,nSrc);
lSig = size(src_sigs,1);
if lMovSig<lSig
    src_mov_sigs = src_sigs(1:lMovSig,:);
else
    for nl=1:ceil(lMovSig/lSig)
        if nl<ceil(lMovSig/lSig)
            src_mov_sigs((nl-1)*lSig +(1:lSig),:) = src_sigs;
        else
            src_mov_sigs((nl-1)*lSig +(1:lMovSig-lSig*(ceil(lMovSig/lSig)-1)),:) = ...
                src_sigs((1:lMovSig-lSig*(ceil(lMovSig/lSig)-1)),:);
        end
    end
end
% get spectra of replicated source signals
winsize = 1024;
% arrival time at RIR sampling points
t_rirs = (0:nRIR-1)*trajectory.rir_dist/trajectory.speed;

shmovsig = 0;
binmovsig = 0;
for ns=1:nSrc
    disp(['Synthesizing motion for source ' num2str(ns)])
    shmovsig = shmovsig + ctf_ltv_direct(src_mov_sigs(:,ns), sh_rirs(:,:,:,ns), t_rirs, scene.fs, winsize);
    
    binmovsig = binmovsig + ctf_ltv_direct(src_mov_sigs(:,ns), bin_rirs(:,:,:,ns), t_rirs, scene.fs, winsize);
end

if WRITE_TO_DISK
    audiowrite([path_out filesep 'moving/orchestra_rev_moving_ref.wav'], shmovsig, scene.fs);
end

end
