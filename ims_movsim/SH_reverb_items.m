clear all
close all

%% Include functions and external tools
gitpath = '../../../git/';
libs = {...
    %'/Higher-Order-Ambisonics/',...
    %'Spherical-Harmonic-Transform/',...
    %'Spherical-Array-Processing/',...
    %'Vector-Base-Amplitude-Panning/'...
    'shoebox-roomsim/'};
for ilib = libs
    addpath(fullfile(gitpath, ilib{1}));
end

addpath('./../')

%render_path = '/work/t405/T40527/chold/rendered/sh_reverb/'
%hrtf_path = '/work/t405/T40527/chold/HRTFs/THK_KU100/'
render_path = './rendered/sh_reverb/'
hrtf_path = '../../../data/HRTFs/THK_KU100/'
mkdir(render_path)
mkdir([render_path '/binaural/']);


%% Parameters
t_60 = 1.2;  % s
fs = 48000;
% SH orders for receivers
rec_orders = [1, 3, 5, 7, 9, 11, 15, 20, 30, 38];

for rec_order = rec_orders

%% Flags
PLOTFLAG = false;
BINFLAG = true;


%% get dry signals
[s_in, fs_in] = audioread('../dry_audio/Clapping-Pink-Guitar-Voice.wav');
assert(fs_in == fs)
%s_in = 0.5*s_in;


%% room definition
room_size = 2*[10.2 7.1 4.2];
fc = [125 250 500 1e3 2e3 4e3 8e3];
nBands = length(fc);
rt60 = t_60 * [1.0 1.0 0.8 0.6 0.45 0.2 0.1];
t = 0:1/fs:2*max(rt60)-1/fs;


%% source positions
% sph(azi, ele, r) in deg
aziEle2aziColat = @(dirs) [dirs(:, 1), pi/2-dirs(:, 2)];

% src_azi = deg2rad([0, 100, -30, -80]);
% src_ele = deg2rad([10, 50, -15, 30]);
%src_r = 0.5*[2, 2.8, 2.2, 1.8];
src_azi = deg2rad([110, -30, 30, 0]);
src_ele = deg2rad([0, 0, 0, 0]);
src_r = [1, 1, 1, 1];


nSrc = length(src_azi);
[src_pos_x, src_pos_y, src_pos_z] = sph2cart(src_azi, src_ele, src_r);
src_pos = [src_pos_x', src_pos_y', src_pos_z'];

% receiver position
rec_pos = 2*[ 4.5 3.3 1.5];
nRec = size(rec_pos,1);



    
%% Apply
sh_ir = generate_sh_reverb(room_size,fc,rt60, src_pos, rec_pos, rec_order, fs);


nSH = (rec_order+1)^2;
nPad = size(sh_ir, 1)-1;
wet_nm = zeros(size(s_in, 1)+nPad, nSH);
    
for iSrc = 1:nSrc
    disp(['Apply reverb to source ' num2str(iSrc)])
%     wet_nm = wet_nm + fftfilt(sh_ir(:, :, iSrc), ...
%                               repmat([s_in(:, iSrc); zeros(nPad, 1)], 1, nSH));
    wet_nm = fftfilt(sh_ir(:, :, iSrc),...
        repmat([s_in(:, iSrc); zeros(nPad, 1)], 1, nSH));

wet_nm = 1/nSrc * wet_nm;

if BINFLAG
    [bin_l, bin_r] = SHtoBin(wet_nm, fs, hrtf_path);
end

    
%% write
outfilename = ['/out', num2str(iSrc), '_rev_SH', num2str(rec_order)];

save([render_path, '/', outfilename, '.mat'], 'wet_nm', 'bin_l', 'bin_r', 'fs', '-v7.3')

audiowrite([render_path, '/', outfilename, 'N3D.wav'], wet_nm, fs, 'BitsPerSample',32);
if BINFLAG
    audiowrite([render_path, '/binaural/', outfilename, 'toBin.wav'], [bin_l, bin_r], fs,'BitsPerSample',32);
end

end

%%
end  % iN


%%
if PLOTFLAG

figure(); plot(t, sh_ir(:, 1, iSrc))


figure()
scatter3(rec(1), rec(2), rec(3), 'k*')
hold on
for i = 1:nSrc
    scatter3(src(i, 1), src(i, 2), src(i, 3))
end
axis equal

set(gca,'XLim',[0 room(1)],'YLim',[0, room(2)],'ZLim',[0 room(3)])
xlabel('x in m'); ylabel('y'); zlabel('z')
legend('rec', 'Clapping', 'Pink', 'Guitar', 'Voice')
hold off;
end