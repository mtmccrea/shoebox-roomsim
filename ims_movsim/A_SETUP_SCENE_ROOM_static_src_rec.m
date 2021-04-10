function [scene, src_sigs] = A_SETUP_SCENE_ROOM_static_src_rec(room, rt60, src_norm, rec_norm, fs, path_out)

% room:     room dimensions in meters [Length Width Height]
% rt60:     rt60 by freq band
% src_norm: specify K sources, in coordinates normalized to the room
%           dimensions [Kx3]
% rec_norm: specify K receivers, in coordinates normalized to the room
%           dimensions [Kx3]

% PATH SETUP
addpath('../')
addpath(genpath('./EXTERNAL_LIBS'))


% GLOBALS
c = 343;

% % LOAD SAMPLES
% [src_sigs, fs_src] = audioread(filepath);
% 
% if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs
% 
% src_sigs = src_sigs./max(abs(src_sigs(:))); 
% 
% src_sigs = src_sigs(:,1:numSources);

% SCENE PARAMETERS
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
% switch roomType
%     case 'anechoic' 
%         room = [8 6 2.5];                   % (L,W,H)
%         rt60 = [0.1 0.1 0.1 0.1 0.1 0.1];   % can't be 0, but the absorption coeffs will be 1 with this
% 
%     case 'medium'  
%         room = [10 7 3];                    % (L,W,H)
%         rt60 = [1.0 0.8 0.7 0.6 0.5 0.4].*0.666;
% 
%     case 'large' 
%         room = [28 22 8];                   % (L,W,H)
%         rt60 = [1.6 1.7 1.6 1.3 1 0.8].*0.666;
% end

%%%%%%%%%%% SOURCE POSITIONS
%room = [35.7 19.8 17.4];
%room = [20 16 7]; %  (L,W,H)

% source positions  
src = src_norm .* room;
% src = src(1:numSources,:); % truncate to requested numSources

% receiver positions (define multiple, if moving)
rec = rec_norm .* room;

% plot scene for visual checking
%plotScene(room, rec, src, 0);

%%%FOR TESTING; DELETE THEM:
%src(1,:) = [rec(1,1) rec(1,2)-2 rec(1,3)]; 
%src(2,:) = [rec(1,1) rec(1,2)+2 rec(1,3)]; 

% reverberation time per octave band
nBands = length(rt60);
% lowest octave band
band_centerfreqs(1) = 125;
for nb=2:nBands, band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end % octave band centerfreqs for RT60
% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room, rt60);
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
scene.AIR_ABS = APPLY_AIR_ABSORPTION;
scene.fs = fs;
scene.c = c;
scene.path_out = path_out;

