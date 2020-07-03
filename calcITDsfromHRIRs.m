function [itd] = calcITDsfromHRIRs(IR, fs, PLOTFLAG)

upsampleIR = upsample(IR, 10);
for i=1:size(upsampleIR,3)
    xcorr_LandR = xcorr(upsampleIR(:,1,i), upsampleIR(:,2,i));
    [~, idx] = max(xcorr_LandR);
    itd(i,1) = (size(upsampleIR,1) -idx)/(fs*10); %#ok
    if itd(i,1) > sqrt(2)/2e3
        itd(i,1) = sqrt(2)/2e3; %#ok
    end
    if itd(i,1) < -sqrt(2)/2e3
        itd(i,1) = -sqrt(2)/2e3; %#ok 
    end
end

if  PLOTFLAG == 1
    itd_ms = itd*1000;
    grid_dirs_deg = IR_pos(1:2,:).';
    grid_dirs = grid_dirs_deg*pi/180;
    grid_dirs2 = [grid_dirs(:,1) pi/2-grid_dirs(:,2)];
    plotSphFunctionTriangle(itd_ms,grid_dirs2,'real') 
end


end

