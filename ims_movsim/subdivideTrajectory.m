function segment_points = subdivideTrajectory(control_points, segment_dist)
%SUBDIVIDETRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

% compute total distance traveled
nCtrlPts = size(control_points,1);
traj_dist = 0;
for npt=1:nCtrlPts-1
    traj_dist = traj_dist + sqrt(sum(abs(control_points(npt+1,:)-control_points(npt,:)).^2,2));
end

% number of rir points
%nSegPts = ceil(traj_dist/segment_dist);
segdist_accum = 0;
segment_points = [];
for npt=1:nCtrlPts-1
    dp_vec = control_points(npt+1,:)-control_points(npt,:);
    dp_uvec = dp_vec/norm(dp_vec);
    dp_tot = sqrt(sum(abs(dp_vec).^2,2));
    pt_new = control_points(npt,:) + segdist_accum*dp_uvec;
    segment_points = [segment_points; pt_new];
    while segdist_accum+segment_dist <= dp_tot
        pt_new = segment_points(end,:) + segment_dist*dp_uvec;
        segdist_accum = segdist_accum+segment_dist;
        segment_points = [segment_points; pt_new];
    end
    if npt<nCtrlPts-1
        segdist_accum = segdist_accum + segment_dist - dp_tot;
    end
end
% add the last point if the one before last is not close enough to the last
% control point
if dp_tot-segdist_accum > segment_dist/2
    segment_points = [segment_points; control_points(end,:)];
end

end