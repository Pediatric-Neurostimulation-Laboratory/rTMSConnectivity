function Feature = ALTO_spTMSEEGConnectivity(EEG, kernel, type, FeMethod, bandname, offset)

% EEG.data = EEG.data(:,:,trial);
% Alto_rsEEGFeatureExtract - Resting EEG feature extraction script.
%                            Features include power spectrum, connectivity
%                            measures (orthogonalized power envelope connectivity, imaginary coherence, wPLI)
%
% Syntax:
%  Feature = Alto_rsEEGFeatureExtract(EEG, kernel, type)
%
% Inputs:
%  EEG         : Clean resting EEG data in EEGLAB structure
%  kernel      : kernel matrix for source localization
%  type        : 'source'  - compute source space features
%                'channel' - compute channel space features
%  NOTE: added by kehli as parameters
%  FeMethod    : 'fourier','powenv','coh','wpli_debiased'
%  bandname    : 'DELTA','THETA','ALPHA','BETA','GAMMA'
%
% Outputs:
%  Feature     : Computed EEG features
%
% Wei Wu, 2019
% Alto Neuroscience, Inc.
%
% Modified by Xiwei She, 2023 to adapt single-pulse TMS-EEG dataset
% Project: rTMS Connectivity

% FeMethod = {'fourier','powenv','coh','wpli_debiased'};
% bandname = {'DELTA','THETA','ALPHA','BETA','GAMMA'};

% NOTE: moved to parameters by kehli
% FeMethod = {'powenv'};
% bandname = {'BETA'};

if strcmp(type,'source')
    sourceflg = 1;
else
    sourceflg = 0;
end
allflg = 0;
cfg.inverseMethod = 'MNE';
cfg.epochSize = 4.0;
all_powspctrm = [];
all_freq = [];
all_coh = [];
all_wpli_debiased = [];


%% computing kernel
if sourceflg == 1
    ds = 3;
    cfg.kernel = kernel;
    sourceData = cfg.kernel*squeeze(mean(EEG.data, 3))*1e-6;
    transmatrix = zeros((size(sourceData,1)/ds), size(sourceData,1));
    for ii = 1: size(sourceData,1)/ds
        sourcetemp = sourceData((ii-1)*ds+1:ii*ds,:);
        R = sourcetemp*sourcetemp';
        [VV, DD] = eig(R);
        [~, I] = sort(diag(DD),'descend');
        transmatrix(ii,(ii-1)*ds+1:ii*ds) = VV(:,I(1))';
    end
    cfg.PCkernel = transmatrix*cfg.kernel;
else
    cfg.PCkernel = kernel;
end
%% Feature calculation
for jj = 1:length(bandname)
    disp(['======' 'Band ' bandname{jj} ' Analysis' '======']);
    for kk = 1:length(FeMethod)
        disp(['======' 'Method ' FeMethod{kk} ' Analysis' '======']);
        
        %% epoch EEG
        cfg.bandname = bandname{jj};
        Data = eeglab2fieldtrip(EEG,'preprocessing','chanlocs');

        % XS: Add offset to remove the artifact of the TMS pulses
        if offset > 0
            currentTimeSpan = [EEG.times(1)/EEG.srate EEG.times(end)/EEG.srate];
            tmpcfg.latency = [offset/EEG.srate, currentTimeSpan(2)];
            Data = ft_selectdata(tmpcfg, Data);
        else
            currentTimeSpan = [EEG.times(1)/EEG.srate EEG.times(end)/EEG.srate];
            tmpcfg.latency = [currentTimeSpan(1), offset/EEG.srate];
            Data = ft_selectdata(tmpcfg, Data);
        end

        switch cfg.bandname
            case 'DELTA'
                cfg.length = 45/1; % length per epoch. Set it in this way because cfg.length*tapsmofrq should be less than 50
            case 'THETA'
                cfg.length = 45/1.5; 
            case 'ALPHA'
                cfg.length = 45/2; 
            case 'BETA'
                cfg.length = 45/8.5; 
            case 'GAMMA'
                cfg.length = 45/9.5; 
        end
        cfg.overlap = 0;
        % Data = ft_redefinetrial(cfg, Data);
        
        %% Feature Calculation Parameters       
        cfg.FeMethod = FeMethod{kk};
        switch cfg.bandname
            case 'DELTA'
                cfg.foi = 1:3; % cfgwpli.tapsmofrq = 1
            case 'THETA'
                cfg.foi = 4:7; % cfgwpli.tapsmofrq = 1.5
            case 'ALPHA'
                cfg.foi = 8:12; % cfgwpli.tapsmofrq = 2
            case 'BETA'
                cfg.foi = 13:30; % cfgwpli.tapsmofrq = 8.5
            case 'GAMMA'
                cfg.foi = 31:50; % cfgwpli.tapsmofrq = 9.5
        end
        
        %% Feature calculation
        switch cfg.FeMethod
            case 'fourier'  % compute band power spectrum
                % set parameters for fourier tranform used in fieldtrip
                cfgfft             = [];
                cfgfft.method      = 'mtmfft';
                cfgfft.taper       = 'dpss';
                cfgfft.output      = 'pow';
                cfgfft.keeptapers  = 'no';
                cfgfft.foi         = cfg.foi;
                cfgfft.channel     = 'all';
                cfgfft.toi         = 0.05:0.05:cfg.epochSize;
                switch cfg.bandname
                    case 'DELTA'
                        cfgfft.tapsmofrq = 1;
                    case 'THETA'
                        cfgfft.tapsmofrq = 1.5;
                    case 'ALPHA'
                        cfgfft.tapsmofrq = 2;
                    case 'BETA'
                        cfgfft.tapsmofrq = 8.5;
                    case 'GAMMA'
                        cfgfft.tapsmofrq = 9.5;
                end
                cfgfft.keeptrials  = 'no';
                cfgfft.pad='nextpow2';      %¼ÓËÙfourier±ä»»
                
                %                 Data = eeglab2fieldtrip(EEG,'preprocessing','chanlocs');
                if sourceflg == 1
                    Data.trial = cellfun(@(s)cfg.PCkernel*s, Data.trial, 'UniformOutput', false);
                    for i = 1:size(cfg.PCkernel,1)
                        Data.label{i} = num2str(i);
                    end
                end
                DataFreqTmp = ft_freqanalysis(cfgfft,Data);
                Feature.(cfg.bandname).powspctrm = DataFreqTmp.powspctrm;
                Feature.(cfg.bandname).freq = cfg.foi;
                clear DataFreqTmp
                all_powspctrm = [all_powspctrm, Feature.(cfg.bandname).powspctrm];
                all_freq = [all_freq, Feature.(cfg.bandname).freq];
                
            case 'powenv'
                cfgpe.bands = 1; % compute powenv at a frequency band
                cfgpe.cpfreq = [cfg.foi(1) cfg.foi(end)];
                cfgpe.bw = 3.8;
                cfgpe.kernel = cfg.PCkernel;
                Feature.(bandname{jj}).powenv = ComputePowerEnvelope(cfgpe,EEG);
                
            case 'coh'
                cfgcoh               = [];
                cfgcoh.method        = 'mtmfft';
                cfgcoh.pad           = 'nextpow2';
                cfgcoh.taper         = 'dpss';
                cfgcoh.output        = 'fourier';     % 'fourier', 'powandcsd'
                cfgcoh.foi           = cfg.foi;
                cfgcoh.channel       = 'all';
                cfgcoh.keeptrials    =  'yes';
                switch cfg.bandname
                    case 'DELTA'
                        cfgcoh.tapsmofrq = 1;
                    case 'THETA'
                        cfgcoh.tapsmofrq = 1.5;
                    case 'ALPHA'
                        cfgcoh.tapsmofrq = 2;
                    case 'BETA'
                        cfgcoh.tapsmofrq = 8.5;
                    case 'GAMMA'
                        cfgcoh.tapsmofrq = 9.5;
                end
                cfgcc.method         = cfg.FeMethod;
                cfgcc.complex        = 'imag';  %'abs' (default), 'angle', 'complex', 'imag', 'real'
                
                % compute features at channel space and then project them to source space
                % compute spectrum
                
                powcsd = ft_freqanalysis(cfgcoh,Data);
                temppw=powcsd.fourierspctrm;
                % source orientation
                for ii=1:size(temppw,3)
                    pwspect(:,:,ii)=temppw(:,:,ii)*cfg.PCkernel';
                end
                powcsd.fourierspctrm=pwspect;
                for lb=1:size(cfg.PCkernel,1) powcsd.label{1,lb}=num2str(lb); end
                clear pwspect temppw;
                tempcoh = ft_connectivityanalysis(cfgcc,powcsd);
                Feature.(bandname{jj}).coh = abs(mean(tempcoh.cohspctrm,3));
                Feature.(bandname{jj}).coh = single(Feature.(bandname{jj}).coh);
                all_coh = cat(3,all_coh,tempcoh.cohspctrm(:,:,1:2:end));
                clear powcsd tempcoh
                
                
            case {'wpli','wpli_debiased'}
                cfgwpli               = [];
                cfgwpli.method        = 'mtmfft';
                cfgwpli.taper         = 'dpss';
                cfgwpli.output        = 'fourier';
                cfgwpli.pad           = 'nextpow2';
                %                     cfgwpli.keeptapers    = 'no';
                cfgwpli.foi           = cfg.foi;
                cfgwpli.channel       = 'all';
                cfgwpli.keeptrials    =  'yes';
                switch cfg.bandname
                    case 'DELTA'
                        cfgwpli.tapsmofrq = 1;
                    case 'THETA'
                        cfgwpli.tapsmofrq = 1.5;
                    case 'ALPHA'
                        cfgwpli.tapsmofrq = 2;
                    case 'BETA'
                        cfgwpli.tapsmofrq = 8.5;
                    case 'GAMMA'
                        cfgwpli.tapsmofrq = 9.5;
                end
                cfgcc.method          = cfg.FeMethod;
                
                cfgcc.complex         = 'complex';  %'abs' (default), 'angle', 'complex', 'imag', 'real'
                
                % compute features at channel space and then project them to source space
                % compute cross spectrum
                
                powcsd = ft_freqanalysis(cfgwpli,Data);
                temppw=powcsd.fourierspctrm;
                for ii=1:size(temppw,3)
                    pwspect(:,:,ii)=temppw(:,:,ii)*cfg.PCkernel';
                end
                powcsd.fourierspctrm=pwspect;
                for i=1:size(cfg.PCkernel,1)
                    powcsd.label{i}=num2str(i);
                end
                clear pwspect temppw;
                tempwpli = ft_connectivityanalysis(cfgcc,powcsd);   % wpli requires >1 trials
                try
                    tempwpli = tempwpli.wplispctrm;
                catch
                    tempwpli = tempwpli.wpli_debiasedspctrm;
                end
                wpli = abs(mean(tempwpli,3));
                %                 wpli = abs(tempwpli);
%                 Feature.(bandname{jj}).wpli_debiased = single(wpli);
                Feature.(bandname{jj}).(cfg.FeMethod) = single(wpli); % Modified by Xiwei
                all_wpli_debiased = cat(3,all_wpli_debiased,abs(tempwpli(:,:,1:2:end)));
                clear powcsd tempwpli wpli
                
                %source orientation
                %             load crospow3003.mat; % for faster computing wpli
                %             sourcepowcsd.label=crospow.label;
                %             sourcepowcsd.labelcmb=crospow.labelcmb;
        end
        
        %         % source ROI-level feature
        %         if sourceflg == 1
        %             if ~strcmp(FeMethod{kk}, 'fourier')
        %                 feature=double(Feature.(bandname{jj}).(FeMethod{kk}));
        %                 feature(isnan(feature)==1) = 0;
        %                 FeatureROI=zeros(length(ROIs),length(ROIs));
        %                 for ii=1:length(ROIs)
        %                     for ll=1:length(ROIs)
        %                         sum=0;
        %                         for iii=1:length(ROIs(ii).Vertices)
        %                             for jjj=1:length(ROIs(ll).Vertices)
        %                                 sum=sum+feature(ROIs(ii).Vertices(iii),ROIs(ll).Vertices(jjj));
        %                             end
        %                         end
        %                         FeatureROI(ii,ll)=sum/(length(ROIs(ii).Vertices)*length(ROIs(ll).Vertices));
        %                         %         wpliROI(j,i)=wpliROI(i,j);
        %                     end
        %                 end
        %                 Feature.ROI.(bandname{jj}).(FeMethod{kk}) = FeatureROI;
        %             end
        %         end
    end
end
if allflg==1
    cfgpe.bands = 1; % compute powenv at a frequency band
    cfgpe.cpfreq = [1 50];
    cfgpe.bw = 3.8;
    cfgpe.kernel = cfg.PCkernel;
    Feature.all.powenv = ComputePowerEnvelope(cfgpe,EEG);
    Feature.all.powspctrm = all_powspctrm;
    Feature.all.freq = all_freq;
    Feature.all.coh = mean(all_coh,3);
    Feature.all.coh = single(Feature.all.coh);
    Feature.all.wpli_debiased = mean(all_wpli_debiased,3);
    Feature.all.wpli_debiased = single(Feature.all.wpli_debiased);
end
end

function [MRO] = ComputePowerEnvelope(peparam,EEG)
for tt=1:size(EEG.data,3)
    if tt==1
        dataConcate=pop_select(EEG,'trial',tt);
    else
        dataTmp=pop_select(EEG,'trial',tt);
        dataConcate=pop_mergeset(dataConcate,dataTmp);
    end
end
EEG=dataConcate;
clear dataConcate;
EEG.data(isnan(EEG.data))=0;

%% Make BPCD
bwInd=0;
clear eegBPCDAll eegBPCD
if peparam.bands==1
    ORIG = EEG;
    EEG = pop_resample(EEG,100); % to determine dimension of eegBPCD
    for f=1:size(peparam.cpfreq,1)
        clear EEG bpcd
        EEG = ORIG;
        EEG.event=[];
        EEG = pop_eegfiltnew(EEG,peparam.cpfreq(f,1),peparam.cpfreq(f,2));
        EEG = pop_resample(EEG,100);
        timepointsLength = size(EEG.data,2);
        eegBPCD = NaN(size(EEG.data,1),timepointsLength,size(peparam.cpfreq,1));
        for c=1:size(EEG.data,1)
            bpcd(c,:) = hilbert(double(EEG.data(c,:)));
        end
        eegBPCD(:,:,f) = bpcd;
    end
    eegBPCDAll(:,:,:)=eegBPCD;
    peparam.BPCD=eegBPCDAll;
    [MRO]=PowerEnvelopeComputation(peparam);
else
    for bw=peparam.bw
        bwInd=bwInd+1;
        ORIG = EEG;
        EEG = pop_resample(EEG,100); % to determine dimension of eegBPCD
        timepointsLength = size(EEG.data,2);
        eegBPCD = NaN(size(EEG.data,1),timepointsLength,length(peparam.cpfreq));
        for f=peparam.cpfreq
            clear EEG bpcd
            if f - bw/2 > 1 && f + bw/2 < 50
                EEG = ORIG;
                EEG = pop_eegfiltnew(EEG,f - bw/2,f + bw/2);
                EEG = pop_reref(EEG,[]);
                EEG = pop_resample(EEG,100);
                for c=1:size(EEG.data,1)
                    bpcd(c,:) = hilbert(double(EEG.data(c,:)));
                end
                eegBPCD(:,:,f) = bpcd;
            end
        end
        eegBPCDAll(:,:,:,bwInd)=eegBPCD;
    end
    clear bpcd;
    peparam.BPCD=eegBPCDAll;
    [MRO]=PowerEnvelopeComputation(peparam);
end
end


function [MRO] = PowerEnvelopeComputation(peparam)
%% Power Envelope - Russ Toll
%%Adapted by Cammie Rolle
cpFreqs=peparam.cpfreq;
kernel=peparam.kernel;
eegBPCD=peparam.BPCD;
bwInd=0;
for cpBandwidth = peparam.bw
    if peparam.bands==1;bwInd=1;else;bwInd=bwInd+1;end
    MRO = single(NaN(size(kernel,1),size(kernel,1),size(cpFreqs,1)));
    for fi=1:size(cpFreqs,1)
        f = cpFreqs(fi);
        eeg = squeeze(eegBPCD(:,:,fi,bwInd));
        [Rplain,Rortho] = orthogonalize(kernel,eeg);
        MRO(:,:,fi,bwInd) = Rortho;
    end
end
end


function [Rplain,Rortho] = orthogonalize(K,BPCD)
% This function calculates the orthogonalized connectivity matrix (Rortho)
% as well as the plain (non-orthogonalized) connectivity matrix (Rplain)
% given a spatial imaging kernel (K) and band-passed complex EEG data (BPCD)
%
% K is a matrix obtained from an inverse solution method such as minimum norm
% estimation using the Brainstorm toolbox. It has dimensionality v vertices x n
% channels.
%
% BPCD is a complex matrix of EEG data that has been band-passed at a certain
% center frequency and bandwidth. For example, EEG data band-passed with center
% frequency 10 Hz and bandwidth 4 Hz to contain the alpha frequency band. This
% data is then made into an analytic signal via the Hilbert transform or using
% Morlet wavelets, etc. It has dimensionality n channels x t time points.

% 1. Input validation.

% 1.1. K and BPCD are required inputs
narginchk(2,2)

% 1.2. K and BPCD must be matrices.
assert(ismatrix(K),'K must be a matrix (vertices x channels).')
assert(ismatrix(BPCD),'BPCD must be a matrix (channels x time points).')

% 1.3. BPCD must be complex-valued.
assert(~isreal(BPCD),'BPCD must be a complex-valued matrix.')

% 1.4. K and BPCD must not have any non-finite values (NaNs).
assert(all(isfinite(K(:))),'Non-finite values detected in K.')
assert(all(isfinite(real(BPCD(:)))),'Non-finite values detected in BPCD.')
assert(all(isfinite(imag(BPCD(:)))),'Non-finite values detected in BPCD.')

% 1.5. Issue a warning if the imaginary component of BPCD is all zero or near
% zero values.
if nansum(nansum(round(imag(BPCD),3))) == 0
    warning('The imaginary component of BPCD appears to contain all zero values, check that your band-passed complex EEG data is correct.')
end

% 2. Calculate source space (SS). Dimensionality is v vertices x t time points.
% Values are computed using single precision as the memory and computation time
% requirements increase substantially when using double precision with
% negligible gain in accuracy.
SS = single(K * BPCD);

% 3. Initialize the directional connectivity matrix (RorthoAB). Dimensionality
% is a square matrix of v vertices x v vertices.
RorthoAB = single(NaN(size(K,1),size(K,1)));

% 4. Determine if a capable GPU is present to perform the orthogonalization.
toolboxInfo = ver;
useGPU = false;
% The Parallel Computing Toolbox is required to use GPU functions.
if any(strcmp({toolboxInfo.Name},'Parallel Computing Toolbox'))
    if gpuDeviceCount > 0
        gpu = gpuDevice(1);
        SSinfo = whos('SS');
        RorthoABinfo = whos('RorthoAB');
        % Check to make sure the detected GPU has enough available memory.
        if gpu.AvailableMemory >= 2 * SSinfo.bytes + RorthoABinfo.bytes
            useGPU = true;
            fprintf('GPU was detected and will be used for calculation.');
        else
            warning('A GPU was detected, but it does not have sufficient memory. Using the CPU to perform the orthogonalization. This will take substantially longer computation time (possibly an hour or more).')
        end
    else
        warning('No GPU was detected, using the CPU to perform the orthogonalization. This will take substantially longer computation time (possibly an hour or more).')
    end
else
    warning('No Parallel Computing Toolbox was detected (required to use GPUs). Using the CPU to perform the orthogonalization. This will take substantially longer computation time (possibly an hour or more).')
end

% 5. Calculate the plain power envelopes (PEplain) by multiplying SS with its
% conjugate. Dimensionality is v vertices x t time points.
PEplain = SS .* conj(SS);

% 6. Calculate the plain connectivity matrix (Rplain) by computing the
% vertex-wise correlations of the natural logarithm of the plain power
% envelopes (PEplain). A small value (tol) is added to the power envelopes to
% prevent the natural logarithm resulting in a non-finite value (i.e., to prevent
% log(0) = -Inf). Dimensionality is a square matrix of v vertices x v vertices.

tol = single(1e-20);
% tol = single(eps);
Rplain = corr(log(PEplain' + tol),log(PEplain' + tol));

% 7. Calculate the directional orthogonalized connectivity matrix (RorthoAB).
% Here, AB refers to the fact that orthogonalization is a directional operation.
% That is, vertex A orthogonalized with respect to vertex B is not the same as
% vertex B orthogonalized with respect to vertex A. This is accomplished by
% orthogonalizing SS with respect to one vertex per iteration, calculating
% the resulting orthogonalized power envelopes (PEortho) and the power envelope
% of the seed vertex (PEseed), and then computing the correlation of the natural
% logarithms of these power envelopes to obtain a measure of connectivity. As
% with the correlation of the natural logarithm of the plain power envelopes, a
% small value (tol) is added to prevent non-finite values.
if useGPU
    SS = gpuArray(SS);
    RorthoAB = gpuArray(RorthoAB);
    tol = gpuArray(tol);
end
% waitbarHandle = waitbar(0,'Initializing','Name','Orthogonalized connectivity');
for iVertex=1:size(SS,1)
    % Seed is the vertex that SS is being orthogonalized with respect to.
    % Dimensionality is 1 x t time points.
    seed = SS(iVertex,:);
    
    % Calculate PEseed by multiplication with its conjugate.
    % Dimensionality is 1 x t time points.
    PEseed = seed .* conj(seed);
    
    % Calculate PEortho by applying Hipp's equation.
    % Dimensionality is v vertices x t time points.
    PEortho = imag(SS .* ( conj(seed) ./ abs(seed) )).^2;
    
    % Calculate the correlation between the natural logarithm of PEseed
    % and the natural logarithm of PEortho. This produces a column vector
    % of v correlation coefficients per iteration.
    
    %    RorthoAB(:,iVertex) = corr(log(PEseed' + tol),log(PEortho' + tol));
    RorthoAB(:,iVertex) = corr(log(PEseed' + tol),log(PEortho' + tol));
    %    RorthoAB(:,iVertex) = corr(PEseed',PEortho');
    
    % Update status bar.
    %    waitbar(iVertex/size(SS,1),waitbarHandle,sprintf('Calculating orthogonalized connectivity.\nVertex %u of %u complete.',iVertex,size(SS,1)));
end
% delete(waitbarHandle);

if useGPU
    RorthoAB = gather(RorthoAB);
    tol = gather(tol);
    reset(gpu);
end

% 8. Calculate the symmetric, corrected orthogonalized connectivity matrix
% (Rortho). For region of interest intersection connectivity analyses, a
% symmetric connectivity matrix is required. Therefore, RorthoAB is averaged
% with its transpose. Due to underestimation of correlation inherent to
% orthogonalization, a correction factor of 0.578499 is applied. See Hipp,
% J.F. et al., 2012. Large-scale cortical correlation structure of spontaneous
% oscillatory activity. Nature Neuroscience, 15(6), pp.884890. for the basis of
% this.
Rortho = ((RorthoAB + RorthoAB') ./ 2) ./ 0.578499;

% 9. Limit the correlation values of Rplain and Rortho to |1 - tol|. This is
% done to prevent non-finite values in the Fisher R to Z transform of the
% correlation coefficients (i.e., fisherz(1) = Inf).
Rplain(Rplain >= 1) = 1 - tol;
Rplain(Rplain <= -1) = -1 + tol;
Rortho(Rortho >= 1) = 1 - tol;
Rortho(Rortho <= -1) = -1 + tol;
end

