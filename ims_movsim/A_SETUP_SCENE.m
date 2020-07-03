function [scene, src_sigs] = A_SETUP_SCENE()

% PATH SETUP
addpath('../')
addpath('stft_conv')
addpath(genpath('./EXTERNAL_LIBS'))
addpath('./resources')
path_out = 'E:/aalto_samples/test_scenes/ORCHESTRA/output';
if ~exist(path_out, 'dir')
    mkdir(path_out)
end

% GLOBALS
fs = 48000;
c = 343;

% LOAD SAMPLES
%[src_sigs, fs_src] = audioread('./samples/bruckner_multichannel.wav');
[src_sigs, fs_src] = audioread('E:/aalto_samples/test_scenes/ORCHESTRA/bruckner_multichannel.wav');
if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs
% remove some silence in the end of the source signals
src_sigs = src_sigs(1:end-fs,:);

%% SCENE PARAMETERS
% Here set up you scene parameters, source locations, microphone locations,
% room boundaries and target reverberation times (or directly wall
% absorption coefficients). If there is no room (anechoic rendering) the
% global origin is arbitrary (e.g. can be at one of the microphones),
% however if there is a room (reverberant rendering), all receiver and
% source positions should be given with respect to the bottom left corner
% of the room (top view), with positive x+ extending to the east, and
% positive y+ extending to the north, while z+ is extending purpendicular
% to them towards the viewer (right-hand rule)
%
%   length/width
%   |----------|
%   ^ y           .
%   |    ^z      /height
%   |   /       /
%   .__/_______.           _
%   | /        |           |
%   |/         |           | width/length
%   o__________.------> x  _
%
% Note that there is no checking for the source-microphone coordinates
% falling inside the boundaries of the room.

% room dimensions
room = [35.7 19.8 17.4]; % some concert hall dims picked up for Musikvereinsaal
% stage height for instrument positions (if needed)
stage_z = 1.5;
% orchestra instrument positions with respect to the center of orchestra
% and height from stage
% 1-14 strings
% 15-18 woodwinds
% 19-24 brass
src2orchestraOrigin_xyz = [
    0        2      1.2
    0        4      1.2
    0        6      1.2
    2        1      1.2
    3        3      1.2
    3        5      1.2
    2       -1      1.2
    3       -3      1.2
    3       -5      1.2
    0       -2      1.2
    0       -4      1.2
    0       -6      1.2
    3       -7      0.4
    5       -6      0.4
    3        1      1.2
    3       -1      1.2
    5        1      1.2
    5       -1      1.2
    5        3      0.5
    5        5      0.5
    6        0      1.2
    6        2      1.2
    6       -2      1.2
    6       -4      0.4];
nSrc = size(src2orchestraOrigin_xyz,1);
% transform to room-based coords for room simulation
orchestraOrigin2roomCorner_xyz = [23 9.5 stage_z];
clear src
src(:,1) = orchestraOrigin2roomCorner_xyz(1) + src2orchestraOrigin_xyz(:,1);
src(:,2) = orchestraOrigin2roomCorner_xyz(2) + src2orchestraOrigin_xyz(:,2);
src(:,3) = orchestraOrigin2roomCorner_xyz(3) + src2orchestraOrigin_xyz(:,3);
% foa/hoa mic positions with respect to array center
d_mic = 4;
mic2arrayOrigin_xyz = [ 0            1.5*d_mic   0;
                        0            0.5*d_mic   0;
                        0           -0.5*d_mic   0;
                        0           -1.5*d_mic   0;
                       -1.25*d_mic   0.5*d_mic   0;
                       -1.25*d_mic  -0.5*d_mic   0];
nMic = size(mic2arrayOrigin_xyz,1);
% array origin with respect to orchestra origin
arrayOrigin2orchestraOrigin_xyz = [-5 0 1.2];
mic2orchestraOrigin_xyz = mic2arrayOrigin_xyz + ones(nMic,1)*arrayOrigin2orchestraOrigin_xyz;
clear rec
rec(:,1) = orchestraOrigin2roomCorner_xyz(1) + mic2orchestraOrigin_xyz(:,1);
rec(:,2) = orchestraOrigin2roomCorner_xyz(2) + mic2orchestraOrigin_xyz(:,2);
rec(:,3) = orchestraOrigin2roomCorner_xyz(3) + mic2orchestraOrigin_xyz(:,3);
% plot scene for visual checking
plotScene(room, rec, src, 0);

% reverberation time per octave band
%rt60 = [2.7706    2.9330    2.9144    2.9332    2.6648    1.8987    1.1478]; % Musikevereinsaal, measured by Lokki and co.
rt60 = 0.8;
nBands = length(rt60);
% lowest octave band
band_centerfreqs(1) = 125;
for nb=2:nBands, band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end % octave band centerfreqs for RT60
% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room, rt60);
% critical distance for the room
[~,  d_critical] = room_stats(room, abs_wall);
APPLY_AIR_ABSORPTION = 1;

% store scene parameters
scene.room = room;
scene.src = src;
scene.rec = rec;
scene.rt60 = rt60;
scene.abs_wall = abs_wall;
scene.abs_bands = band_centerfreqs;
scene.AIR_ABS = APPLY_AIR_ABSORPTION;
scene.fs = fs;
scene.c = c;
scene.path_out = path_out;

