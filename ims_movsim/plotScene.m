function [h_fig, h_mic, h_src, h_tri] = plotScene(room, rec_xyz, src_xyz, TOP_VIEW_ONLY)
%PLOTSCENE Summary of this function goes here
%   Detailed explanation goes here

if nargin<4, TOP_VIEW_ONLY = 1; end

h_fig = figure;
nMic = size(rec_xyz,1);
nSrc = size(src_xyz,1);

% plot room top view
if ~TOP_VIEW_ONLY, subplot(211), end
if all(room~=0)
    rectangle('Position',[0,0,room(1),room(2)],...
         'LineWidth',2,'LineStyle','-')
end

% plot sources & receivers
h_src = line(src_xyz(:,1),src_xyz(:,2),'linestyle','none','marker','*','markersize',14,'color','b','linewidth',3);
h_mic = line(rec_xyz(:,1),rec_xyz(:,2),'linestyle','none','marker','o','markersize',14,'color','r','linewidth',3);
for ns=1:nSrc, text(src_xyz(ns,1),src_xyz(ns,2),['  s' num2str(ns)],'fontsize',14), end
for nm=1:nMic, text(rec_xyz(nm,1),rec_xyz(nm,2),['  m' num2str(nm)],'fontsize',14), end

% plot triangulation
if nMic==2
    h_tri = line(rec_xyz(:,1),rec_xyz(:,2),'linestyle','--','color','k');
elseif nMic>2
    tri = delaunay(rec_xyz(:,1), rec_xyz(:,2));
    hold on, h_tri = triplot(tri, rec_xyz(:,1), rec_xyz(:,2),'linestyle','--','color','k');
elseif nMic==1 % quick fix for single-mic case
    h_tri = line([rec_xyz(:,1), rec_xyz(:,1)],[rec_xyz(:,2),rec_xyz(:,2)], ...
        'linestyle','--','color','k');
end
grid, axis equal
xlabel('x (m)'), ylabel('y (m)')

if ~TOP_VIEW_ONLY
    subplot(212)
    if all(room~=0)
        rectangle('Position',[0,0,room(1),room(3)],...
            'LineWidth',2,'LineStyle','-')
    end
    % plot sources & receivers
    line(src_xyz(:,1),src_xyz(:,3),'linestyle','none','marker','*','markersize',14,'color','b','linewidth',3)
    line(rec_xyz(:,1),rec_xyz(:,3),'linestyle','none','marker','o','markersize',14,'color','r','linewidth',3)
    for ns=1:nSrc, text(src_xyz(ns,1),src_xyz(ns,3),['  s' num2str(ns)],'fontsize',14), end
    for nm=1:nMic, text(rec_xyz(nm,1),rec_xyz(nm,3),['  m' num2str(nm)],'fontsize',14), end
    grid, axis equal
    xlabel('x (m)'), ylabel('z (m)')
end   

end
