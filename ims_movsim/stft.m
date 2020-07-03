function spectrum = stft(insig, winsize, fftsize)

lSig = size(insig,1);
nCHin = size(insig,2);

% time-frequency processing
if nargin<3 || isempty(fftsize)
    fftsize = 2*winsize;
end
hopsize = winsize/2;
nBins = fftsize/2 + 1;
nWindows = ceil(lSig/winsize);
nFrames = 2*nWindows+1;

% zero pad the signal's start and end for STFT
insig_pad = [zeros(winsize/2, nCHin); insig; zeros(nFrames*winsize/2-lSig, nCHin)];
clear insig

spectrum = zeros(nBins, nFrames, nCHin);
% transform window (hanning)
x = 0:(winsize-1);
win = sin(x.*(pi/winsize))'.^2;

% processing loop
idx = 1;
nf = 1;

while nf <= nFrames
    % Window input and transform to frequency domain
    insig_win = win*ones(1,nCHin) .* insig_pad(idx+(0:winsize-1),:);
    inspec = fft(insig_win, fftsize);
    inspec = inspec(1:nBins,:); % keep up to nyquist
    spectrum(:,nf,:) = inspec;
    
    % advance sample pointer
    idx = idx + hopsize;
    nf = nf + 1;
end

end
