function [binsig, closest_dirs_deg] = binauralize(insig, sig_dirs_deg, hrirs, hrir_dirs_deg, fs_in, fs_hrir )
%BINAURALIZE Summary of this function goes here
%   Detailed explanation goes here

nDirs = size(sig_dirs_deg,1);
nCH = size(insig,2);
lSig = size(insig,1);
if (nCH~=nDirs), error('signal DoAs and number of channels do not match'); end

[idx_closest, closest_dirs_rad] = findClosestHRIRs(hrir_dirs_deg*pi/180, sig_dirs_deg*pi/180);
%idx_closest = knnsearch( unitSph2cart(hrir_dirs_deg*pi/180), unitSph2cart(sig_dirs_deg*pi/180) );
closest_dirs_deg = closest_dirs_rad*180/pi;

if (fs_in ~= fs_hrir)
    disp(['Resampling to ' num2str(fs_hrir) ' to match HRIR samplerate'])
    insig = resample(insig, fs_hrir, fs_in);
end

hrirs_closest = hrirs(:,:,idx_closest);
binsig = zeros(lSig, 2);
for n=1:2
    binsig(:,1) = sum( fftfilt(squeeze(hrirs_closest(:,1,:)), insig), 2);
    binsig(:,2) = sum( fftfilt(squeeze(hrirs_closest(:,2,:)), insig), 2);
end

% normalize
maxa = max(max(abs(binsig)));
binsig = 0.5*binsig/maxa;

end

function [idx_closest, dirs_closest, angle_diff] = findClosestHRIRs(hrir_dirs_rad, src_dirs_rad)

nHRIR = size(hrir_dirs_rad,1);
nDirs = size(src_dirs_rad,1);

xyz_hrir = unitSph2cart(hrir_dirs_rad);
xyz_src = unitSph2cart(src_dirs_rad);

idx_closest = zeros(nDirs,1);
for nd=1:nDirs
    [~, idx_closest(nd)] = max(dot(xyz_hrir, repmat(xyz_src(nd,:),nHRIR,1), 2));
end
    
dirs_closest = hrir_dirs_rad(idx_closest,:);
angle_diff = acos( dot(xyz_hrir(idx_closest,:), xyz_src, 2) );

end
