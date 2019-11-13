classdef FluorRecon
    
    % Multicolor F19 Reconstruction Code
    % Usage: R=FluorRecon
    
    % AMSTERDAM UMC, 2019
    %
    % JASPER SCHOORMANS
    % GUSTAV STRIJKERS
    %
    %
    
    
    
    properties
        P                           % Parameters
        Recon                       % Reconstructions
        Recon_originalFOV
        LinearRecon                 % Linear reconstructions
        Functions                   % Functions
        Data                        % Raw k-space data
        k                           % Data vector
        app                         % App with graphical user interface
    end
    
    
    
    methods
        
        function R = FluorRecon(varargin)
            
            R = R.DefineP;
            
            if nargin == 1
                R.P.scanfolders=varargin{1};
            
            else %no folders given, use GUI    %%%%% Gustav: OBSOLETE FOR THE APP %%%%%
            
                ii=1;
                while  ii<5
                    fld = uigetdir('Choose base folder of scan - (up to 4)');
                    if fld~=[0]
                        R.P.scanfolders{ii}=fld;
                        ii=ii+1;
                    else
                        break;
                    end
                    
                end
            end
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=loadBrukerFiles(R)
            
            % Load the data for the FLASH_MP_19F sequence
            % Fully sampled data
            
            TextMessage(R.app,"Loading scan ... ");
            TextMessage(R.app,"Single file scan ... ");
            
            [K,ScanParams] = ReadFIDSequential(R.P.scanfolders{1},R.P.Seq.ndirs,R.P.directions);
            
            R.Data.K = K;
            R.P.ScanParams = ScanParams;
            
        end
        
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=loadBrukerFilesCS(R)
            
            % Load the data for the FLASH_MP_CS_19F sequence
            % Undersampled data
            
            TextMessage(R.app,"Loading undersampled scan ... ");
            TextMessage(R.app,"Single file scan ... ");
            
            % kspace filling pattern
            R.app.TextMessage('User-defined k-space filling ... ');
            dimy = R.app.brukerinfo2.acq.size(2);
            dimz = R.app.brukerinfo2.acq.size(3);
            ky1rep = round(R.app.brukerinfo2.acq.spatial_phase_1 * (dimy/2-0.5) + dimy/2);
            kz1rep = round(R.app.brukerinfo2.acq.spatial_phase_2 * (dimz/2-0.5) + dimz/2);
            
            % For now assuming no_repetitions = 1
            no_repetitions = 1;
            ky=[]; for i=1:no_repetitions ky = [ky, ky1rep]; end
            kz=[]; for i=1:no_repetitions kz = [kz, kz1rep]; end
            
            [K,ScanParams] = ReadFIDSequentialCS(R.P.scanfolders{1},R.P.Seq.ndirs,R.P.directions,ky,kz);
            
            R.Data.K = K;
            R.P.ScanParams = ScanParams;
            
        end
        
        
        
        % -----------------------------------------------------------------------------
        
        
        
        
        function R = ConjugateGradientRecon(R)
            
            TextMessage(R.app,"CG recon ... ");
            
            R = toGPU(R);
            nouter = 3;  % number of outer iterations
            
            % outer iterations
            for outeriter = 1:nouter
                
                TextMessage(R.app,strcat("Outer iteration ",num2str(outeriter)," ... "));
                
                if outeriter == 1
                    RCG=zeros([2*R.P.N,1]);
                end
                
                % inner iterations
                RCG = nl_conjgrad_fluor_3D(R.app,R.Functions.M, R.k, RCG,R.P.niter,nouter,outeriter,R.Functions.TVOP2,R.P.lambda,[R.P.nx R.P.ny R.P.nz],R.P.visualization,R.P.visualization_slice);
                
            end
            
            R.Recon{2,1}=gather(R.Functions.rr(abs(RCG(1:R.P.N))));
            R.Recon{2,2}=gather(R.Functions.rr(abs(RCG(R.P.N+1:2*R.P.N))));
            
            TextMessage(R.app,"FINISHED ... ");
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=BlindDeconvRecon(R)
            
            TextMessage(R.app,"Blind Deconvolution ... ");
            
            M = R.Functions.M;
            nouter = 3;
            
            for outeriter = 1:nouter
                
                TextMessage(R.app,strcat("Outer iteration ",num2str(outeriter)," ... "));
                
                if outeriter==1 
                    RCG=zeros([2*R.P.N,1]); 
                end
                
                RCG = nl_conjgrad_fluor_3D(R.app,R.Functions.M, R.k, RCG,R.P.niter,nouter,outeriter,R.Functions.TVOP2,R.P.lambda,[R.P.nx R.P.ny R.P.nz],R.P.visualization,R.P.visualization_slice);
                
                kk = BlindDeconvolutionOptimization(RCG,R.P.N,R.LinearRecon,R.P.ny,R.P.Spectrum2,R.Functions.rr,R.P.directions,1,R.P.visualization);
                
                UpdatedSpectrum=cell(1,4);
                UpdatedSpectrum{1}=kk(1:R.P.nx).';
                UpdatedSpectrum{2}=kk(R.P.nx+1:2*R.P.nx).';%k(nx+1:2*nx);
                UpdatedSpectrum{3}=kk(2*R.P.nx+1:3*R.P.nx).';%k(2*nx+1:3*nx);
                UpdatedSpectrum{4}=kk(3*R.P.nx+1:4*R.P.nx).';%k(3*nx+1:4*nx);
                
                M= MakeMeasurementOperator(UpdatedSpectrum,R.P.nx,R.P.ny,R.P.nz,R.P.directions);
                
            end
            
            R.Recon{2,1}=R.Functions.rr(abs(RCG(1:R.P.N)));
            R.Recon{2,2}=R.Functions.rr(abs(RCG(R.P.N+1:2*R.P.N)));
            
            TextMessage(R.app,"FINISHED ... ");
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=PseudoInverse(R)
            
            TextMessage(R.app,"Pseudo-inverse Reconstruction ...");
            
            first_guess=pinv(R.Functions.M)*R.k;
            
            R.Recon{1,1}=R.Functions.rr(abs(first_guess(1:R.P.N)));
            R.Recon{1,2}=R.Functions.rr(abs(first_guess(R.P.N+1:2*R.P.N)));
            
            TextMessage(R.app,"FINISHED ... ");
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=ApplyCheckBoard(R)
            
            if R.P.CheckBoard
                
                TextMessage(R.app,"Applying checkerboard pattern ... ");
                
                for ii=1:R.P.nacq
                    ch=create_checkerboard(size(R.Data.K{ii}));
                    R.Data.K{ii}=bsxfun(@times,R.Data.K{ii},ch);
                end
                
            end
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = VectorizeData(R)
            
            cmd='R.k=[';
            
            for ii = 1:R.P.nacq
                
                cmd = [cmd,'R.Functions.vec(R.Data.K{',num2str(ii),'})'];
                
                if ii < R.P.nacq
                    cmd = [cmd,'; '];
                else
                    cmd = [cmd,'];'];
                end
                
            end
            
            eval(cmd);
        
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = PhaseRemoval(R)
            
            if R.P.phaseremoval
                
                TextMessage(R.app,"Removing phase of image data ... ");
                
                R.k = R.Functions.FB*abs(opInverse(R.Functions.FB)*R.k); 
            
            end
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = LinRecon(R)
            
            TextMessage(R.app,"Linear reconstruction ... ");
            
            R.LinearRecon=pinv(R.Functions.FB)*R.k;    %linear recon of data
            
            TextMessage(R.app,"FINISHED ... ");
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = FixDoubleFOV(R)
            
            TextMessage(R.app,"Fixing double FOV ... ");
            
            % dealing with double FOV in readout direction (2 options)
            if ~R.P.doubleFOV
                
                TextMessage(R.app,"Removing odd k-space points ... ");
                
                for i=1:R.P.nacq
                    
                    if R.P.ScanParams{i}.AntiAlias(1)==2;
                        
                        Ktemp=R.Data.K{i};
                    
                        if strcmp(R.P.directions{i},'FH') | strcmp(R.P.directions{i},'HF')
                            Ktemp=Ktemp(:,1:2:end,:); %begin at 1 or 2? wrt reversing direction...
                            
                        elseif strcmp(R.P.directions{i},'RL') | strcmp(R.P.directions{i},'LR')
                            Ktemp=Ktemp(1:2:end,:,:);
                        else
                            TextMessage(R.app,"ERROR: Direction not implemented ... ");
                        end
                        
                        R.Data.K{i}=Ktemp;
                        
                    end
                    
                end
                
            else
                
                for i=1:R.P.nacq
                    
                    TextMessage(R.app,"Doubling FOV in two phase encoding dimensions... ");
                    
                    fft3_meas = @(I,dim) fftshift(fftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
                    ifft3_meas = @(I,dim) fftshift(ifftn(ifftshift(I))); % iFFT in measurement direction  %%? which dimension??
                    
                    Ktemp=R.Data.K{i};
                
                    if strcmp(R.P.directions{i},'FH') | strcmp(R.P.directions{i},'HF')
                        
                        pad = (size(Ktemp,2)-size(Ktemp,1))/2;
                        I1=ifft3_meas(Ktemp);
                        I1=padarray(I1,[pad 0 pad]);
                        Ktemp=fft3_meas(I1);
                        
                    elseif strcmp(R.P.directions{i},'RL') | strcmp(R.P.directions{i},'LR')
                        
                        pad=(size(Ktemp,1)-size(Ktemp,2))/2;
                        I1=ifft3_meas(Ktemp);
                        I1=padarray(I1,[0 pad pad]);
                        Ktemp=fft3_meas(I1);
                        
                    else
                        
                        TextMessage(R.app,"ERROR: Direction not implemented ... ");
                    
                    end
                    
                    R.Data.K{i}=Ktemp;
                                        
                    % BW/pix essentially doubles
                
                end
                
            end
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = DefineP(R)
            % Defines standard settings for Fluor Recon
            
            R.P.doubleFOV = 0;                                  % no double FOV for now
            R.P.dim = 3;                                        % 3D data
            R.P.phaseremoval = 1;                               % option to remove phase
            R.P.phasecheckerboardremoval = ~R.P.phaseremoval;
            R.P.CheckBoard = 1;                                 % apply checkerboard to k-space to remove readout/phase encoding phases
            R.P.use_empirical_alpha = 1;                        % overwrite theoretical peaks with empirical
            R.P.TranslationCorr = 1;                            % automatically translate k-spaces using image registration
            % R.P.ZeroFill = 0;                                 % determined by the app. enlarge kspace (resolution) by ZF
            
            R.P.visualization = 0;                              % visualization only in the App
            
            % R.P.BW                                            % Acquisition bandwidth, defined in the app
            % R.P.ppm                                           % parts per million / defines the field strength
            
            R.Data.K=[];                                        % raw k-space
            R.Data.improton = [];                               % proton image
            
            R.Recon_originalFOV = cell(3,2);
            R.Recon=cell(3,2);
            
            R.P.Brukerfolder = '';
            R.P.visualization_slice = 1;
            
            % R.P.Seq.Sequential                                % Sequential measurement (multiple dirs in one file)
            % R.P.Seq.ndirs                                     % number of readout directions
            
            % R.P.useGPU=0;                                     % using GPU (faster)
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = UpdateP(R)
            % Updates parameters and operators after data loading
            TextMessage(R.app,"Update operators and derive parameters from data ...");
            
            R.P.nacq=size(R.Data.K,2);                          % number of separate acquisitions
            assert(R.P.nacq==size(R.P.directions,2),'Directions parameter not ok');
            
            R.P.phasecheckerboardremoval=~R.P.phaseremoval;
            
            % size dependent parameters
            R.P.nx=R.P.ScanParams{1}.Enc(1)/R.P.ScanParams{1}.AntiAlias(1); % imsizes
            R.P.ny=R.P.ScanParams{1}.Enc(2);                    % imsizes
            R.P.nz=R.P.ScanParams{1}.Enc(3);                    %imsizes
            
            if R.P.doubleFOV
                R.P.nx=R.P.nx*2;
                R.P.ny=R.P.ny*2;
                R.P.nz=R.P.nz*2;
            end
            
            if R.P.ZeroFill
                
                for ii=1:R.P.nacq
                    [nx,ny,nz]=size(R.Data.K{ii}); %before removing oversampling - so recalc
                    R.Data.K{ii}=padarray(R.Data.K{ii},[nx/2,ny/2,nz/2]);
                end
                
                R.P.nx=2*R.P.nx; R.P.ny=2*R.P.ny; R.P.nz=2*R.P.nz;
            
            end
            
            TextMessage(R.app,"TO DO: automate number of pixel setting (oversampling) ... ");
            
            R.P.N=R.P.nx*R.P.ny*R.P.nz; % number of pixels in one image
            
            %after fix double FOV
            R.Recon{1,1}=zeros(R.P.nx,R.P.ny,R.P.nz);
            R.Recon{1,2}=zeros(R.P.nx,R.P.ny,R.P.nz);
            R.Recon{2,1}=zeros(R.P.nx,R.P.ny,R.P.nz);
            R.Recon{2,2}=zeros(R.P.nx,R.P.ny,R.P.nz);
            R.Recon{3,1}=zeros(R.P.nx,R.P.ny,R.P.nz);
            R.Recon{3,2}=zeros(R.P.nx,R.P.ny,R.P.nz);
            
            if R.P.doubleFOV
                R.P.BWpix=R.P.BW/(R.P.nx/2); % for a double FOV, the BW/pix is essentially doubled
            else
                R.P.BWpix=R.P.BW/(R.P.nx);
            end
            
            %functions
            R.Functions.vec= @(I) reshape(I,[numel(I), 1]); %vectorize
            R.Functions.rr2 = @(I) reshape(I,[R.P.nx,R.P.ny,R.P.nz,2]); %resize to two images (C2/C1)
            R.Functions.rr1 = @(I) reshape(I,[R.P.nx,R.P.ny,R.P.nz,R.P.nacq]); %resize to nr of acqs.
            R.Functions.rr = @(I) reshape(I,[R.P.nx,R.P.ny,R.P.nz]);  %resize to 1 image
            
            % operatorsFSp
            R.Functions.F3=opDFT3(R.P.nx,R.P.ny,R.P.nz,1);
            
            TextMessage(R.app,"Constructing TV operator ... ");
            R.Functions.TVOP2=MakeTVOperator(R.P.nx);
            
            % block diagonal operator
            inp = [];
            for ii = 1:R.P.nacq
                inp = [inp,'R.Functions.F3,'];
            end
            R.Functions.FB=eval(['opBlockDiag(',inp(1:end-1),')']);  % 4 fourier ops - one for each image
            
            % add Fourier operator for all nacq (only for changing when
            % undersampling/different resolutions)
            for ii=1:R.P.nacq
                R.Functions.F{ii} = R.Functions.F3;
            end
           
            % 24-10-2019
            % TO DO: ADD HERE 
            % 
            
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = CalculateSpectra(R)
            
            % Gustav: Adapted. Now includes peak positions and heights from 2 files, which are read in the app
            
            TextMessage(R.app,"Calculating the spectrum of compound 1 and compound 2 ... ");
            
            % Compound 1 peak positions and peak heights
            C1 = (R.P.ppm.*R.P.empirical_ppm1)./R.P.BWpix;
            C1_alpha = R.P.empirical_alpha1;
            
            pixlocs1=round(1-C1);
            pixlocs1(pixlocs1<0)=R.P.nx+pixlocs1(pixlocs1<0);
            R.P.Spectrum1=zeros(1,R.P.nx);
            R.P.Spectrum1(pixlocs1)=C1_alpha./sum(C1_alpha(:));
            
            % Compound 2 peak positions and peak heights
            C2 = (R.P.ppm.*R.P.empirical_ppm2)./R.P.BWpix;
            C2_alpha = R.P.empirical_alpha2;
            
            pixlocs2=round(1-C2);
            pixlocs2(pixlocs2<0)=R.P.nx+pixlocs2(pixlocs2<0);
            R.P.Spectrum2=zeros(1,R.P.nx);
            R.P.Spectrum2(pixlocs2)=C2_alpha./sum(C2_alpha(:));
       
            % construct the measurement operator
            for ii=1:R.P.nacq
                R.P.SpectrumC1{ii}=R.P.Spectrum1;
            end
            
            for ii=1:R.P.nacq
                R.P.SpectrumC2{ii}=R.P.Spectrum2;
            end
            
            R.Functions.M = MakeMeasurementOperator2C(R.P.SpectrumC1,R.P.SpectrumC2,R.P.nx,R.P.ny,R.P.nz,R.P.directions,R.Functions.F);
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = Undersample(R)
            % Retrospectively undersample k-space
            % Experimental 
            
            TextMessage(R.app,"K-space undersampling (compressed sensing) ... ");
            
            pattern=bart('poisson -Y64 -Z64 -C20 -y1.2 -z1.2 ');
            undersampling=1./(sum(pattern(:))/numel(pattern));
            
            % apply pattern (for all directions)
            K=R.Functions.rr1(R.k);
            datatemp=[] ;
            pattern1=repmat(pattern,[64 1 1]); %FH/HF
            pattern2=repmat(permute(pattern,[2 1 3]),[1 64 1]); %LR/RL
            
            for i = 1:R.P.nacq
                ktemp = K(:,:,:,i);
                if strcmp(R.P.directions{i},'FH') | strcmp(R.P.directions{i},'HF')
                    ktemp=ktemp(pattern1(:)==1);
                elseif strcmp(R.P.directions{i},'RL') | strcmp(R.P.directions{i},'LR')
                    ktemp=ktemp(pattern2(:)==1);
                else
                    error('Direction not implemented')
                end
                datatemp = [datatemp;ktemp(:)]; %save and reshape back to vector
            end
            
            R.k=datatemp;
            
            %%%update Fourier Operators to undersampled Fourier operators
            p=[];
            for i = 1:R.P.nacq
                if strcmp(R.P.directions{i},'FH') | strcmp(R.P.directions{i},'HF')
                    newF{i} = opExcise(R.Functions.F3,~pattern1(:),'rows');
                elseif strcmp(R.P.directions{i},'RL') | strcmp(R.P.directions{i},'LR')
                    newF{i} = opExcise(R.Functions.F3,~pattern2(:),'rows');
                end
            end
            
            %update block diagonal Fourier Operator
            inp = [];
            for ii = 1:R.P.nacq
                inp=[inp,'newF{',num2str(ii),'},'];
            end
            R.Functions.FB = eval(['opBlockDiag(',inp(1:end-1),')']);  % 4 fourier ops - one for each image
            R.Functions.F = newF;

        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R = TestMeasOp(R)
            % Test constructed Measurement operator on dummy data
            % Purpose: to see if peak heights/directions correspond to
            % linear recon of acquired image
            
            TextMessage(R.app,"Testing dummy data ... ");
            
            %used to compare the measurement operator with real data
            temp=zeros(R.P.nx,R.P.nx,R.P.nx,2);
            temp(R.P.nx/2,R.P.nx/2 -1,R.P.nx/2,1)=5;
            temp(R.P.nx/2,R.P.nx/2 -1,R.P.nx/2,1)=5;
            temp(R.P.nx/2,R.P.nx/2 +1,R.P.nx/2,1)=5;
            temp(R.P.nx/2 +1,R.P.nx/2,R.P.nx/2,1)=5;
            temp(R.P.nx/2 +1,R.P.nx/2 +1,R.P.nx/2,1)=5;
            
            temp2=R.Functions.M*temp(:);
            temp3=pinv(R.Functions.FB)*temp2;
            temp4=(R.Functions.rr1(temp3));
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=TranslationCorr(R)
            % Correct for offcenter acq. by automatic registration and translation in ksp
            
            if R.P.TranslationCorr
                
                TextMessage(R.app,"Autocorrection for off-center acquisition ... ");
                
                for ii =2:R.P.nacq
                    
                    TextMessage(R.app,strcat("Registering direction ",num2str(ii)," ... "));
                    
                    R.Data.K{ii} = registration_correction_19F(R.app,R.Data.K{ii},R.Data.K{1},0);
                end
            end
            
        end
        
        
        % -----------------------------------------------------------------------------
        
        
        function R=toGPU(R)
            % Moves arrays to the GPU for faster computation 
            
            if R.P.useGPU
                
                TextMessage(R.app,"Moving k-space to GPU ... ");
                
                R.k=gpuArray(R.k);
                
            end
        end
        
        
        % -----------------------------------------------------------------------------
        
        
    end
end
