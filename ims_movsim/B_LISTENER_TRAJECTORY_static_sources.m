function trajectory = B_LISTENER_TRAJECTORY_1source_static(scene, GENERATE_ANIM)

if nargin<2, GENERATE_ANIM=0; end
rec = scene.rec;
nMic = size(rec,1);

%% Moving receiver parameters
% Here the trajectory followed by the listener is defined. We assume
% a constant moving speed for the simulations, and with the listener moving
% in line segments on the x-y plane of the microphone array.
%
% In the following example the listener is moving in a path connecting
% directly microphones (on the edges of the triangulation connecting the
% mics)

% define trajectory control points
cPts(1,:) = rec(1,:);  
trajectory.cPts = cPts;
% moving speed
mov_speed = 0.5; % m/s
% compute total distance traveled
nCPts = size(cPts,1);
traj_dist = 0;
for npt=1:nCPts-1
    traj_dist = traj_dist + sqrt(sum(abs(cPts(npt+1,:)-cPts(npt,:)).^2,2));
end
% traveling time
mov_time = traj_dist/mov_speed;

trajectory.speed = mov_speed;
trajectory.dist = traj_dist;
trajectory.time = mov_time;

%% Generate animation
% Here an animation video of the motion is created if desired, in order to
% display during listening tests, as a visual aid
if GENERATE_ANIM
    % generate all intermediate points for the animation
    framerate = 30; %fps
    frametime = 1/framerate;
    framedist = mov_speed*frametime;
    p_vid = cPts;%subdivideTrajectory(cPts, framedist);
    
    % create video
    [~,h_mic] = plotScene(scene.room, scene.rec, scene.src);
    x_head = [2 cos((10:350)*pi/180) 2];
    y_head = [0 sin((10:350)*pi/180) 0];
    g = hgtransform;
    patch('XData',0.5*x_head,'YData',0.5*y_head,'FaceColor','m','Parent',g)
    set(gca, 'xlim', [-1 scene.room(1)+1], 'ylim', [-1 scene.room(2)+1])
    set(gcf, 'Position', [65 292 907 647])
    markercolor = h_mic.Color;
    markersize = h_mic.MarkerSize;
    v = VideoWriter([scene.path_out filesep 'animation.avi'], 'MPEG-4');
    open(v);
    for np = 1:size(p_vid,1)
        p_current = p_vid(np,:);
        [~, closestmic] = min(sum(abs(ones(nMic,1)*p_current - rec).^2,2));
        rec_closest = rec(closestmic, :);
        h_lin = line(rec_closest(1), rec_closest(2), ...
            'linestyle','none','marker','o','markersize', markersize, ...
            'MarkerEdgeColor', 'none', 'MarkerFaceColor', markercolor);
        
        g.Matrix = makehgtform('translate', p_current);
        title(['time: ' num2str(round((np-1)*frametime)) ' sec'])
        drawnow
        frame = getframe(gcf);
        writeVideo(v,frame);
        delete(h_lin)
    end
    close(v);
end

%% Trajectory RIR points
% Define points in the trajectory for computation of static room impulse 
% responses for the approximation of the time-variant reverberation as the
% listener moves.

% distance betwen successive RIR computation points on the trajectory
trajectory.rir_dist = 0.2;
% number of rir points
trajectory.rir_pts = cPts;%subdivideTrajectory(cPts, trajectory.rir_dist);

%% Trajectory points for parameter interpolation
% Here we compute the points and the interpolation weights for each
% processing frame of the synthesized rendering using the surrounding 
% microphones

% number of samples for computation of parameters
framesize = 2048;
frametime = framesize/scene.fs;
framedist = frametime*mov_speed;
trajectory.intrp_pts = cPts;%subdivideTrajectory(cPts, framedist);
% compute interpolation weights for every segment point in the trajectory
trajectory.intrp_weights = 1;%getPlanarTriInterpWeights(scene.rec, trajectory.intrp_pts);

end