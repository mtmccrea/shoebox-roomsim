function [scene, src_sigs] = A_SETUP_SCENE(roomType, filepath, numSources, path_out)

% PATH SETUP
addpath('../')
addpath('stft_conv')
addpath('fromChris/')
addpath(genpath('./EXTERNAL_LIBS'))
%addpath('./resources')


% GLOBALS
fs = 48000;
c = 343;

% LOAD SAMPLES
%[src_sigs, fs_src] = audioread('./samples/bruckner_multichannel.wav');
%[src_sigs, fs_src] = audioread('/Users/mccorml1/Documents/aalto_samples/test_scenes/ORCHESTRA/bruckner_multichannel.wav');
%[src_sigs, fs_src] = audioread('/Users/mccorml1/Documents/aalto_samples/test_scenes/4SPEAKERS_maleEng_femEng_maleDan_femDan.wav');
[src_sigs, fs_src] = audioread(filepath);

if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs

src_sigs = src_sigs./max(abs(src_sigs(:))); 
src_sigs = src_sigs(:,1:numSources);

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

%%%%%%%%%%% ROOM 
switch roomType
    case 'anechoic' 
        room = [8 6 2.5];
        rt60 = [0.1 0.1 0.1 0.1 0.1 0.1]; % can't be 0, but the absorption coeffs will be 1 with this

    case 'medium'  
        room = [10 7 3];
        rt60 = [1.0 0.8 0.7 0.6 0.5 0.4].*0.666;

    case 'large' 
        room = [28 22 8];
        rt60 = [1.6 1.7 1.6 1.3 1 0.8].*0.666;
end

%%%%%%%%%%% SOURCE POSITIONS
%room = [35.7 19.8 17.4];
%room = [20 16 7]; %  (L,W,H)
% source positions  
clear src 
src(1,:) = [room(1)*0.8 room(2)*0.7 room(3)*0.35]; 
if numSources >=2
    src(2,:) = [room(1)*0.3 room(2)*0.2 room(3)*0.45]; 
end
if numSources >=3
    src(3,:) = [room(1)*0.15 room(2)*0.8 room(3)*0.25]; 
end
if numSources >=4
    src(4,:) = [room(1)*0.7 room(2)*0.3 room(3)*0.85];  
end

% receiver positions (define multiple, if moving)
clear rec 
rec = [room(1)*0.1 room(2)*0.45 room(3)*0.4; ...
       room(1)*0.5 room(2)*0.4 room(3)*0.2; ...
       room(1)*0.85 room(2)*0.5 room(3)*0.55; ...
       room(1)*0.6 room(2)*0.7 room(3)*0.3; ...
       room(1)*0.35 room(2)*0.7 room(3)*0.65;];
% plot scene for visual checking
plotScene2(room, rec, src, 0);

% reverberation time per octave band 
nBands = length(rt60);
abs_wall_ratios = [0.75 0.86 0.56 0.95 0.88 1];
% lowest octave band
band_centerfreqs(1) = 125;
for nb=2:nBands, band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end % octave band centerfreqs for RT60
% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room, rt60);%,abs_wall_ratios);
% critical distance for the room
%[~,  d_critical] = room_stats(room, abs_wall);
APPLY_AIR_ABSORPTION = 1;

% store scene parameters
scene.room = room;
scene.src = src;
scene.rec = rec;
scene.rt60 = rt60;
scene.abs_wall = abs_wall;
scene.abs_bands = band_centerfreqs;
%scene.abs_wall_ratios = abs_wall_ratios;
scene.AIR_ABS = APPLY_AIR_ABSORPTION;
scene.fs = fs;
scene.c = c;
scene.path_out = path_out;




