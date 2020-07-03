function convsig = ctf_ltv_direct(sig, irs, ir_times, fs, winsize)

% STFT
hopsize = winsize/2;
fftsize = 2*winsize;
nBins = fftsize/2+1;

% IRs
lIr = size(irs,1);
if length(size(irs))==2
    nIrs = size(irs,2);
    nCHir = 1;
elseif length(size(irs))==3
    nIrs = size(irs,3);
    nCHir = size(irs,2);
end
if nIrs ~= length(ir_times)
    error('asd')
end
% number of STFT frames for the IRs (half-window hopsize)
nIrWindows = ceil(lIr/winsize);
nIrFrames = 2*nIrWindows+1;
% number of STFT frames for the signal (half-window hopsize)
lSig = size(sig,1);
nSigWindows = ceil(lSig/winsize);
nSigFrames = 2*nSigWindows+1;
nCHsig = size(sig,2);

% quantize the timestamps of each IR to multiples of STFT frames (hopsizes)
tStamps = round((ir_times*fs+hopsize)/hopsize);
% create the two linear interpolator tracks, for the pairs of IRs between
% timestamps
nIntFrames = tStamps(end);
Gint = zeros(nIntFrames, nIrs);
for ni=1:nIrs-1
    tpts = tStamps(ni):tStamps(ni+1);
    ntpts = length(tpts);
    Gint(tpts,ni) = 1-(0:ntpts-1)/(ntpts-1);
    Gint(tpts,ni+1) = (0:ntpts-1)/(ntpts-1);
end

% compute spectra of irs
if nCHir ==1, irspec = zeros(nBins, nIrFrames, nIrs);
else, irspec = zeros(nBins, nIrFrames, nCHir, nIrs);
end
for ni=1:nIrs
    if nCHir ==1
        irspec(:,:,ni) = stft(irs(:,ni), winsize);
    else
        irspec(:,:,:,ni) = stft(irs(:,:,ni), winsize);
    end
end

% compute input signal spectra
sigspec = stft(sig, winsize);

% initialize interpolated time-variant ctf
Gbuf = zeros(nIrFrames,nIrs);
if nCHir == 1,  ctf_ltv = zeros(nBins,nIrFrames);
else,           ctf_ltv = zeros(nBins,nIrFrames,nCHir);
end
S = zeros(nBins, nIrFrames);

nFrames = min(nSigFrames, nIntFrames);
% processing loop
idx = 1;
nf = 1;

convsig = zeros(winsize/2 + nFrames*winsize/2 + fftsize-winsize, nCHir);
%inspec_pad = [inspec zeros(nBins,nIrFrames)];
inspec_pad = sigspec;
while nf <= nFrames
    % compute interpolated ctf
    Gbuf(2:end,:) = Gbuf(1:end-1,:);
    Gbuf(1,:) = Gint(nf,:);
    if nCHir ==1
        for nif = 1:nIrFrames
            ctf_ltv(:,nif) = squeeze(irspec(:,nif,:))*Gbuf(nif,:).';
        end
    else
        for nch=1:nCHir
            for nif = 1:nIrFrames
                ctf_ltv(:,nif,nch) = squeeze(irspec(:,nif,nch,:))*Gbuf(nif,:).';
            end
        end
    end
    
    inspec_nf = inspec_pad(:,nf);
    S(:,2:nIrFrames) = S(:,1:nIrFrames-1);
    S(:,1) = inspec_nf;
    
    convspec_nf = squeeze(sum(repmat(S,[1 1 nCHir]).*ctf_ltv,2));
    convspec_nf = [convspec_nf; conj(convspec_nf(end-1:-1:2,:))];
    convsig_nf = ifft(convspec_nf,fftsize,1);

    % overlap-add synthesis
    convsig(idx+(0:fftsize-1),:) = convsig(idx+(0:fftsize-1),:) + convsig_nf;
    % advance sample pointer
    idx = idx + hopsize;
    nf = nf + 1;
    disp([num2str(round(nf/nFrames*100)) '%'])
end
convsig = convsig(winsize+1:nFrames*winsize/2,:);

end