classdef T4ResolveBuilder0 < mlfourdfp.IT4ResolveBuilder
	%% T4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Mar-2016 21:29:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	
    
    properties (Constant)
        FRAC_FRAME_THRESH = 0.2
    end
    
    properties
        doMaskBoundaries = false
        firstCrop = 0.5
        keepForensics = true
        logger % TODO
        msktgenThresh = 0
        NRevisions 
        
        frameRegLog     
        t4ResolveLog
    end
    
	properties (Dependent)
        cFrameIndices
        fortranFrameIndices
        frame0
        frameF
        frames
        targetFrame
        
        blurArg
        buildVisitor        
        gaussArg
        inverseCrop
        petBlur
 		product
        referenceImage
        referenceWeight
        resolveTag
        sessionData
        sourceImage
        sourceWeight
        studyData
        xfm
    end
    
    methods %% GET/SET
        function g    = get.blurArg(this)
            if (~isempty(this.blurArg_))
                g = this.blurArg_;
                return
            end
            
            %g = 5.5; % Avi's empirical suggestion for MR, PET alignment
            g = this.petBlur * 1.2790697674; 
        end
        function this = set.blurArg(this, s)
            assert(isnumeric(s));
            this.blurArg_ = s;
        end
        function this = set.buildVisitor(this, v)
            assert(isa(v, 'mlfourdfp.FourdfpVisitor'));
            this.buildVisitor_ = v;
        end
        function v    = get.buildVisitor(this)
            v = this.buildVisitor_;
        end
        function c    = get.cFrameIndices(this)
            c = this.fortranFrameIndices;
            c = c - 1;
        end
        function f    = get.fortranFrameIndices(this)
            indices = 1:length(this.frames);
            f = indices(this.frames ~= 0);
        end
        function g    = get.frame0(this)
            for fr = 1:length(this.frames)
                if (this.frames(fr))
                    g = fr;
                    return
                end
            end
        end
        function g    = get.frameF(this)
            for fr = length(this.frames):-1:1
                if (this.frames(fr))
                    g = fr;
                    return
                end
            end
        end
        function g    = get.frames(this)
            if (~isempty(this.frames_))
                g = this.frames_;
                return
            end
            g = ones(1, this.readLength(this.sessionData.tracerNative('typ', 'fqfp')));
        end
        function this = set.frames(this, s)
            assert(isnumeric(s));
            maxFrames = this.readLength(this.sessionData.tracerNative('typ', 'fqfp'));
            if (isempty(s))
                this.frames_ = [];
                return
            end
            if (length(s) > maxFrames)
                this.frames_ = s(1:maxFrames);
            end
            if (length(s) <= maxFrames)
                s_ = zeros(1, maxFrames);
                s_(1:length(s)) = s;
                this.frames_ = s_;
            end
        end
        function g    = get.inverseCrop(this)
            inv = round(1/this.firstCrop);
            g = [inv inv 1];
        end
        function g    = get.gaussArg(this)
            %g = 1.1;
            g = this.F_HALF_x_FWHM/this.blurArg;
        end
        function pb   = get.petBlur(this)
            pb = this.sessionData.petBlur;
        end
        function prod = get.product(this)
            prod = this.product_;
        end
        function rw   = get.referenceWeight(~)
            rw = [];
        end
        function ri   = get.referenceImage(~)
            ri = [];
        end
        function g    = get.resolveTag(this)
            g = this.resolveTag_;
        end
        function this = set.resolveTag(this, s)
            assert(ischar(s));
            this.resolveTag_ = s;
            this.sessionData_.t4ResolveBuilder = this;
        end
        function sd   = get.sessionData(this)
            sd = this.sessionData_;
        end
        function si   = get.sourceImage(this)
            si = this.sessionData.tracerNative;
        end
        function sw   = get.sourceWeight(~)
            sw = [];
        end
        function g    = get.studyData(this)
            g = this.sessionData.studyData;
        end
        function g    = get.targetFrame(this)
            g = this.targetFrame_;
        end
        function this = set.targetFrame(this, s)
            assert(isnumeric(s));
            if (1 <= s && s <= length(this.frames))
                this.targetFrame_ = s;
            end
        end
        function x    = get.xfm(this)
            x = this.t4Path;
        end
    end
    
    methods (Static)
        function tf     = completed(sessd)
            assert(isa(sessd, 'mlpipeline.SessionData'));
            this = mlfourdfp.T4ResolveBuilder0('sessionData', sessd);
            tf = lexist(this.completedTouchFile, 'file');
        end
        function fp     = coulombPotential(fp, varargin)
            import mlfourdfp.*;
            ip = inputParser;
            addRequired( ip, 'fp', @(x) lexist([fp '.4dfp.ifh']));
            addParameter(ip, 'precision', T4ResolveBuilder0.coulombPrecision, @isnumeric);
            parse(ip, fp, varargin{:});
            
            ic  = mlfourd.ImagingContext([fp '.4dfp.ifh']);
            ic  = ic.coulombPotential('precision', ip.Results.precision);
            fp  = ic.fileprefix;
            icc = ic.clone;
            icc.saveas([fp '.4dfp.ifh']);
        end
        function pr     = coulombPrecision(varargin)
            ip = inputParser;
            addParameter(ip, 'precision',  0.125,           @isnumeric);
            addParameter(ip, 'fileprefix', '',              @ischar);
            addParameter(ip, 'mmppix',     [1 1 1]*2.08626, @isnumeric);
            addParameter(ip, 'blur',       [],              @isnumeric);
            parse(ip, varargin{:});
            
            if (~isempty(ip.Results.blur))
                mmppix = mean(abs(ip.Results.mmppix));
                if (~isempty(ip.Results.fileprefix))
                    mmppix = mean(abs( ...
                             mlfourdfp.T4ResolveBuilder0.readScalingFactors(ip.Results.fileprefix)));
                end
                pr = (0.6875 / ip.Results.blur) * (mmppix / 2.08626); % empirically similar to blur of 5.5 mm
                return
            end
            pr = ip.Results.precision;
        end
        function          diaryv(name)
            assert(ischar(name));
            if (~isempty(getenv('PRINTV')))
                diary(sprintf('mlfourdfp_T4ResolveBuilder_%s.txt', name)); 
            end
        end
        function fn     = fourdfpImg(fp)
            fn = [fp '.4dfp.img'];
        end
        function fp     = frameFileprefix(fp, fr)
            fp = sprintf('%s_frame%i', fp, fr);
        end
        function f      = frameNumber(str, offset)
            names = regexp(str, '\w+(-|_)(F|f)rame(?<f>\d+)', 'names');
            f = str2double(names.f) + offset;
        end
        function          printv(varargin)
            if (~isempty(getenv('PRINTV')))
                fprintf('mlfourdfp.T4ResolveBuilder0.');
                fprintf(varargin{:});
                fprintf('\n');
            end
        end
        function mmppix = readScalingFactors(varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp',  @ischar);
            parse(ip, varargin{:});
            
            mmppix = zeros(1,3);
            for mu = 1:3
                [~,f] = mlbash(sprintf('awk ''/scaling factor \\(mm\\/pixel\\) \\[%i\\]/{print $NF}'' %s.4dfp.ifh', mu, ip.Results.fqfp));
                mmppix(mu) = str2double(f);
            end
        end
        function f      = readLength(varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp',  @ischar);
            addOptional(ip, 'mu', 4, @isnumeric); % tensor index, or "dimension"
            parse(ip, varargin{:});
            [~,f] = mlbash(sprintf('awk ''/matrix size \\[%i\\]/{print $NF}'' %s.4dfp.ifh', ip.Results.mu, ip.Results.fqfp));
            f = str2double(f);
        end
    end
    
	methods
 		function this = T4ResolveBuilder0(varargin)
 			%% T4RESOLVEBUILDER
 			%  @param named buildVisitor
            %  @param named sessionData
            %  @param named frames
            %  @param named nonEmptyFrames
            %  @param named targetFrame
            %  @param named NRevisions
            %  @param named resolveTag
 			
            this = this@mlfourdfp.IT4ResolveBuilder;
            if (0 == nargin); return; end
            
            %% invoke copy-ctor
            
            if (1 == nargin && isa(varargin{1}, 'mlfourdfp.T4ResolveBuilder0'))
                aCopy = varargin{1};
                
                % properties
                this.NRevisions = aCopy.NRevisions;
                this.doMaskBoundaries = aCopy.doMaskBoundaries;
                this.firstCrop = aCopy.firstCrop;
                this.frames = aCopy.frames;
                this.framesToSkip = aCopy.framesToSkip;
                this.inverseCrop = aCopy.inverseCrop;
                this.keepForensics = aCopy.keepForensics;
                this.msktgenThresh = aCopy.msktgenThresh;
                this.resolveTag = aCopy.resolveTag;
                this.targetFrame = aCopy.targetFrame;
                
                % properties (Access = protected)
                this.blurArg_ = aCopy.blurArg_;
                this.buildVisitor_ = aCopy.buildVisitor_;
                this.frameRegLog = aCopy.frameRegLog;
                this.frameRegOrdinalFrameF_ = aCopy.frameRegOrdinalFrameF_;
                this.frameRegOrdinalFrame0_ = aCopy.frameRegOrdinalFrame0_;
                this.mprage_ = aCopy.mprage_;
                this.nonEmptyFrames_ = aCopy.nonEmptyFrames_;
                this.product_ = aCopy.product_;
                this.referenceImage_ = aCopy.referenceImage_;
                this.referenceWeight_ = aCopy.referenceWeight_;
                this.sessionData_ = aCopy.sessionData_;
                this.t4ResolveLog = aCopy.t4ResolveLog;
                this.trash_ = aCopy.trash_;
                this.xfm_ = aCopy.xfm_; 
                return
            end
            
            %% manage parameters 
            
            import mlfourdfp.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'buildVisitor',   FourdfpVisitor, @(x) isa(x, 'mlfourdfp.FourdfpVisitor'));
            addParameter(ip, 'sessionData',    [],             @(x) isa(x, 'mlpipeline.ISessionData'));
            addParameter(ip, 'frames',         [],             @isnumeric);
            addParameter(ip, 'nonEmptyFrames', [],             @isnumeric);
            addParameter(ip, 'targetFrame',    1,              @isnumeric);
            addParameter(ip, 'NRevisions',     1,              @isnumeric);
            addParameter(ip, 'resolveTag',     'resolved',     @ischar);
            parse(ip, varargin{:});
            
            this.buildVisitor_ = ip.Results.buildVisitor;
            this.sessionData_ = ip.Results.sessionData;
            this.sessionData_.t4ResolveBuilder = this;
            if (~isempty(ip.Results.nonEmptyFrames))
                this.nonEmptyFrames_ = ip.Results.nonEmptyFrames;
                this.frames = this.nonEmptyFrames_;
            end
            if (~isempty(ip.Results.frames))
                this.frames = ip.Results.frames;
            end
            this.targetFrame = ip.Results.targetFrame;
            this.NRevisions = ip.Results.NRevisions;
            this.resolveTag = ip.Results.resolveTag;
        end
        
        function this = resolve(this, varargin)
            %% RESOLVE iteratively calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param destMask   "
            %  @param source     "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param useMetricGradient:  cf. ${TRANSFER}/cross_modal_intro.pdf
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            import mlfourdfp.*;
            ip = inputParser;
            addParameter(ip, 'dest',       '',              @ischar);
            addParameter(ip, 'source',     '',              @FourdfpVisitor.lexist_4dfp);
            addParameter(ip, 'destMask',   'none',          @ischar);
            addParameter(ip, 'sourceMask', 'none',          @ischar);
            addParameter(ip, 'framesMask', 'maskOrdFrames', @ischar);
            addParameter(ip, 'destBlur',   this.blurArg,    @isnumeric); % fwhh/mm
            addParameter(ip, 'sourceBlur', this.blurArg,    @isnumeric); % fwhh/mm            
            addParameter(ip, 'NRevisions', this.NRevisions, @isnumeric);
            addParameter(ip, 'firstCrop',  this.firstCrop,  @isnumeric); % [0 1]
            addParameter(ip, 'frames',     this.frames,     @isnumeric); % 1 to keep; 0 to skip
            addParameter(ip, 'log',        '/dev/null',     @ischar);
            addParameter(ip, 'rnumber',    this.sessionData.rnumber, @isnumeric);
            addParameter(ip, 't40',        this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'atlas',      this.atlas('typ', 'fp'), @ischar);
            addParameter(ip, 'par',        false,           @islogical);
            addParameter(ip, 'single',     false,           @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;       
            ipr = this.expandFrames(ipr); 
            ipr = this.findFrameBounds(ipr);
            ipr = this.expandBlurs(ipr);
            ipr.resolved = ipr.source; % initialize this.iterate
            while (ipr.rnumber <= ipr.NRevisions)
                ipr = this.iterate(ipr, ...
                    ipr.resolved, ...
                    this.iterationFileprefix(ipr.dest, ipr.rnumber));
                assert(ipr.rnumber < 10);
            end            
            
            this.buildVisitor.imgblur_4dfp(ipr.resolved, this.blurArg);
            this.product_ = ipr.resolved;
            this.teardownResolve(ipr);
        end
        function ipr  = iterate(this, ipr, source, dest)
            ipr.source = source;
            ipr.dest = dest;            
            ipr = this.cropOrCopy(ipr); % crop/copy ipr.source to ipr.dest
            this.frameRegLog = loggerFilename( ...
                dest, 'func', 'T4ResolveBuilder_frameReg', 'path', this.logPath);
            this.t4ResolveLog = loggerFilename( ...
                dest, 'func', 'T4ResolveBuilder_t4ResolveAndPaste', 'path', this.logPath);
            
            ipr = this.frameReg(ipr);
            ipr = this.t4ResolveAndPaste(ipr); 
            ipr = this.teardownIterate(ipr);
            ipr.rnumber = ipr.rnumber + 1;
        end
        function ipr  = frameReg(this, ipr)
            if (ipr.par)
                ipr = this.frameRegPar(ipr);
                return
            end
            if (ipr.single)
                ipr = this.frameRegSingle(ipr);
                return
            end
            ipr = this.frameRegSerial(ipr);
        end
        function ipr  = t4ResolveAndPaste(this, ipr)
            %% T4RESOLVEANDPASTE
            %  @param ipr is a struct containing fields dest, frames, rnumber
            dest = mybasename(ipr.dest);
            imgFns = this.frameFileprefix(dest, this.targetFrame);
            for f = 1:length(ipr.frames)
                if (ipr.frames(f))
                    if (f ~= this.targetFrame)
                        imgFns = [imgFns ' ' this.frameFileprefix(dest, f)]; %#ok<AGROW>
                    end
                end
            end                     
            
            this.buildVisitor.t4_resolve(this.resolveTag, imgFns, ...
                'options', '-v -m -s', ...
                'log', this.t4ResolveLog);
            this.t4imgFrames(ipr, this.resolveTag);
            this.pasteFrames(ipr, this.resolveTag);

            ipr.resolved = sprintf('%s_%s', dest, this.resolveTag); 
            movefile([this.resolveTag '.mat0'], [ipr.resolved '_' datestr(now, 30) '.mat0']);
            movefile([this.resolveTag '.sub'],  [ipr.resolved '_' datestr(now, 30) '.sub']);
        end
        function        teardownLogs(this)
            if (this.keepForensics); return; end
            
            ensuredir(this.logPath);
            movefiles('*.log', this.logPath); 
            movefiles('*.txt', this.logPath);   
            movefiles('*.lst', this.logPath);    
            movefiles('*.mat0', this.logPath);   
            movefiles('*.sub', this.logPath); 
        end
        function        teardownT4s(this)
            if (this.keepForensics); return; end
            
            ensuredir(this.t4Path);
            movefiles('*_t4', this.t4Path);
            sessd = this.sessionData;
            movefile( ...
                fullfile(this.t4Path, ...
                    sprintf('%s_to_%s_t4', sessd.mpr('typ', 'fp'), sessd.atlas('typ', 'fp'))), ...
                pwd);
        end
        function ipr  = teardownIterate(this, ipr)
            this.teardownLogs;
            this.teardownT4s;
        end
        function        teardownResolve(this, ipr)
            if (this.keepForensics); return; end
            
            for r = 1:this.NRevisions
                
                fp0 = this.iterationFileprefix(ipr.dest, r);
                for f = 1:length(ipr.frames)
                    if (ipr.frames(f))
                        %delete(sprintf('%s_frame%i.4dfp.*', fp0, f));
                        delete(sprintf('%s_frame%i_b*.4dfp.*', fp0, f));
                        delete(sprintf('%s_frame%i_C*.4dfp.*', fp0, f));
                        %delete(sprintf('%s_frame%i_%s.4dfp.*', fp0, f, this.resolveTag));
                        delete(sprintf('%s_frame%i_g*.nii.gz', fp0, f));
                    end
                end

                sessd = this.sessionData;
                sessd.rnumber = r;
                %delete(sprintf('%s_frame*.4dfp.*', sessd.tracerRevision('typ','fqfp')));
            end
            
            delete(sprintf('%s_*_*.4dfp.*', ipr.framesMask));
        end
        function pth  = logPath(this)
            pth = fullfile(this.sessionData.tracerLocation, 'Log', '');
            if (~isdir(pth))
                mkdir(pth);
            end
        end
        function pth  = onAtlasPath(this)
            pth = fullfile(this.sessionData.tracerLocation, 'Atlas', '');
            if (~isdir(pth))
                mkdir(pth);
            end
        end
        function pth  = t4Path(this)
            pth = fullfile(this.sessionData.tracerLocation, 'T4', '');
            if (~isdir(pth))
                mkdir(pth);
            end
        end
        
        %% UTILITY 
        
        function a         = atlas(this, varargin) 
            a = this.sessionData.atlas(varargin{:});
        end
        function sessd     = adjustedSessionData(this, varargin)
 			%% ADJUSTEDSESSIONDATA
 			%  @param [param-name, param-value[, ...]]
            %         'ac'          is logical
            %         'rnumber'     is numeric
            %         'sessionPath' is a path to the session data
            %         'snumber'     is numeric
            %         'tag'         is appended to the fileprefix
            %         'tracer'      is char
            %         'vnumber'     is numeric

            ip = inputParser;
            sessd = this.sessionData;
            addParameter(ip, 'ac',          sessd.attenuationCorrected, @islogical);
            addParameter(ip, 'rnumber',     sessd.rnumber,              @isnumeric);
            addParameter(ip, 'sessionPath', sessd.sessionPath,          @isdir);
            addParameter(ip, 'snumber',     sessd.snumber,              @isnumeric);
            addParameter(ip, 'tag',         sessd.tag,                  @ischar);
            addParameter(ip, 'tracer',      sessd.tracer,               @ischar);
            addParameter(ip, 'vnumber',     sessd.vnumber,              @isnumeric);
            parse(ip, varargin{:});
            
            sessd.attenuationCorrected = ip.Results.ac;
            sessd.rnumber              = ip.Results.rnumber;
            sessd.sessionPath          = ip.Results.sessionPath;
            sessd.snumber              = ip.Results.snumber;
            sessd.tag                  = ip.Results.tag;
            sessd.tracer               = ip.Results.tracer;
            sessd.vnumber              = ip.Results.vnumber;
        end
        function ipr       = findFrameBounds(~, ipr)
            ipr.frame0 = 1;
            ipr.frameF = length(ipr.frames);
            for f = 1:length(ipr.frames)
                if (ipr.frames(f))
                    ipr.frame0 = f;
                    break
                end
            end
            for f = length(ipr.frames):-1:ipr.frame0
                if (ipr.frames(f))
                    ipr.frameF = f;
                    break
                end
            end
        end
        function fqfn      = completedTouchFile(this)
            fqfn = fullfile(this.sessionData.tracerNACLocation, ...
                sprintf('.%s_%s_completed.touch', ...
                    lower(this.sessionData.tracer), class(this)));
        end
        function ipr       = cropOrCopy(this, ipr)
            if (1 == ipr.rnumber)
                try
                    if (~this.buildVisitor.lexist_4dfp(ipr.dest))
                        this.buildVisitor.cropfrac_4dfp(ipr.firstCrop, ipr.source, ipr.dest);
                    end
                catch ME
                    handexcept(ME);
                end
                return
            end
            this.buildVisitor.copy_4dfp( ...
                this.resolvedFileprefix( ipr.dest, ipr.rnumber-1), ...
                this.iterationFileprefix(ipr.dest, ipr.rnumber));
            this.copyResolvedToNewRevision(ipr);
        end
        function             copyResolvedToNewRevision(this, ipr)
            %% COPYRESOLVEDTONEWREVISION opportunistically reuses existing frame files from the last iteration
            
            dt = mlsystem.DirTool( ...
                sprintf('%s_frame*_%s.4dfp.ifh', this.sessionData.tracerRevision('typ','fp'), this.resolveTag));
            for f = 1:length(dt.fns)
                resolvedFrame = mybasename(dt.fns{f});
                idx_   = regexp(resolvedFrame, '_frame\d+');
                idx__  = regexp(resolvedFrame, sprintf('_%s$', this.resolveTag));
                newRev = [ipr.dest resolvedFrame(idx_:idx__-1)];
                try
                    this.buildVisitor.move_4dfp(resolvedFrame, newRev);
                catch ME
                    handexcept(ME);
                end
            end
        end
        function [s,r]     = ensureNifti(this, varargin)
            %% ENSURENIFTI
            %  @param filename is any string descriptor found in an existing file on the filesystem;
            %  ensureNifti will search for files with extensions .nii, .nii.gz or .4dfp.ifh.
            %  @returns s, the bash status; r, any bash messages.  ensureNifti ensures files are .nii.gz.
            
            ip = inputParser;
            addRequired(ip, 'filename', @ischar);
            parse(ip, varargin{:});
            
            s = 0; r = '';
            if (2 == exist(ip.Results.filename, 'file'))
                if (lstrfind(ip.Results.filename, '.nii'))
                    return
                end
                if (lstrfind(ip.Results.filename, '.mgz'))
                    fp = myfileprefix(ip.Results.filename);
                    [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii.gz', fp, fp));
                    return
                end
                [s,r] = this.buildVisitor.nifti_4dfp_ng(myfileprefix(ip.Results.filename));
                assert(lexist(myfilename(ip.Results.filename), 'file'));
                return
            end
            if (2 == exist([ip.Results.filename '.nii'], 'file'))
                [s,r] = this.ensureNifti([ip.Results.filename '.nii']);
                return
            end
            if (2 == exist([ip.Results.filename '.nii.gz'], 'file'))
                return
            end
            if (2 == exist([ip.Results.filename '.mgz'], 'file'))
                [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii.gz', ip.Results.filename, ip.Results.filename));
                return
            end    
            if (2 == exist([ip.Results.filename '.4dfp.ifh'], 'file'))
                if (2 == exist([ip.Results.filename '.nii'], 'file'))
                    return
                end
                [s,r] = this.ensureNifti([ip.Results.filename '.4dfp.ifh']);
                return
            end
            error('mlfourdfp:fileNotFound', ...
                  'T4ResolveBuilder0.ensureNifti could not find files of form %s', ip.Results.filename);            
        end
        function [s,r]     = ensure4dfp(this, varargin)
            %% ENSURE4DFP
            %  @param filename is any string descriptor found in an existing file on the filesystem;
            %  ensureNifti will search for files with extensions .4dfp.ifh.
            %  @returns s, the bash status; r, any bash messages.  ensure4dfp ensures files are .4dfp.ifh.            
            
            ip = inputParser;
            addRequired(ip, 'filename', @ischar);
            parse(ip, varargin{:});
            
            s = 0; r = '';
            if (2 == exist(ip.Results.filename, 'file'))
                if (lstrfind(ip.Results.filename, '.4dfp'))
                    return
                end
                if (lstrfind(ip.Results.filename, '.mgz'))
                    fp = myfileprefix(ip.Results.filename);
                    [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii', fp, fp)); %#ok<ASGLU>
                    [s,r] = this.ensure4dfp([fp '.nii']);
                    return
                end
                [s,r] = this.buildVisitor.nifti_4dfp_4(myfileprefix(ip.Results.filename));
                assert(lexist(myfilename(ip.Results.filename, '.4dfp.ifh'), 'file'));
                return
            end
            if (2 == exist([ip.Results.filename '.4dfp.ifh'], 'file'))
                return
            end
            if (2 == exist([ip.Results.filename '.nii'], 'file'))
                [s,r] = this.ensure4dfp([ip.Results.filename '.nii']);
                return
            end
            if (2 == exist([ip.Results.filename '.nii.gz'], 'file'))
                [s,r] = this.ensure4dfp([ip.Results.filename '.nii.gz']);
                return
            end      
            if (2 == exist([ip.Results.filename '.mgz'], 'file'))
                [s,r] = mlbash(sprintf('mri_convert %s.mgz %s.nii', ip.Results.filename, ip.Results.filename)); %#ok<ASGLU>
                [s,r] = this.ensure4dfp([ip.Results.filename '.nii']);
                return
            end            
            error('mlfourdfp:fileNotFound', ...
                  'T4ResolveBuilder0.ensureNifti could not find files of form %s', ip.Results.filename);     
        end
        function ipr       = expandBlurs(~, ipr)
            assert(isstruct(ipr));
            ipr.sourceBlur = ipr.sourceBlur .* ipr.frames;
            ipr.destBlur   = ipr.destBlur   .* ipr.frames;
        end
        function ipr       = expandFrames(this, ipr)
            assert(isstruct(ipr));
            if (length(ipr.frames) == 1)
                lenSrc = this.readLength(ipr.source);
                if (ipr.frames > 0)
                    ipr.frames = ones(1, lenSrc);
                    return
                end
                ipr.frames = zeros(1, lenSrc);
            end
        end
        function fqfp      = fileprefixBlurred(this, fqfp)
            fqfp = this.buildVisitor.fileprefixBlurred(fqfp, this.blurArg);
        end
        function fqfp      = fileprefixGaussed(this, fqfp)
            fqfp = this.buildVisitor.fileprefixGaussed(fqfp, this.blurArg);
        end
        function fqfp      = fileprefixMskt(~, fqfp)
            fqfp = [fqfp '_mskt'];
        end
        function fqfp      = fileprefixMsktgen(this, fqfp)
            fqfp = this.fileprefixMskt( ...
                this.fileprefixGaussed(this.fileprefixSumt(fqfp)));
        end
        function fqfp      = fileprefixSumt(this, fqfp)
            assert(ischar(fqfp));
            fqfp = sprintf('%s_sumt%i', fqfp, sum(this.frames));
        end
        function fqfp      = fileprefixTmp(this, varargin)
            ip = inputParser;
            addParameter(ip, 'path', this.sessionData.sessionPath, @isdir);
            parse(ip, varargin{:});
            
            fqfp = fullfile(ip.Results.path, sprintf('tmp_%s', datestr(now, 30)));            
        end
        function t4c       = framesToResolvedT4s(this, fp, findices)
            assert(ischar(fp));
            assert(isnumeric(findices));
            fp = mybasename(fp);
            t4c = cellfun( ...
                @(x) fullfile(this.t4Path, sprintf('%s_frame%i_to_%s_t4', fp, x, this.resolveTag)), ...
                num2cell(findices), ...
                'UniformOutput', false);
        end
        function fqfps     = lazyBlurredFrames(this, extractedFps, ipr)
            ip = inputParser;
            addRequired(ip, 'extractedFps', @iscell);
            addRequired(ip, 'ipr', @isstruct);
            parse(ip, extractedFps, ipr);
            
            fqfps = extractedFps;
            f___ = 0;
            for f = 1:length(ipr.frames)                
                if (ipr.frames(f) && all(ipr.sourceBlur(f) > 0))
                    f___ = f___ + 1;
                    fqfps{f___} = this.blurredFrame(extractedFps{f___}, 'blur', ipr.sourceBlur(f));
                end
            end
        end
        function fqfp      = lazyExtractFrame(this, varargin)
            ip = inputParser;
            addRequired(ip, 'dest', @mlfourdfp.FourdfpVisitor.lexist_4dfp);
            addRequired(ip, 'fr', @isnumeric);
            addOptional(ip, 'fqfp', '', @ischar);
            parse(ip, varargin{:});
            if (isempty(ip.Results.fqfp))
                fqfp = this.frameFileprefix(ip.Results.dest, ip.Results.fr);
            else
                fqfp = ip.Results.fqfp;
            end            
            this.buildVisitor.extract_frame_4dfp(ip.Results.dest, ip.Results.fr, ['-o' fqfp]);
        end
        function fqfps     = lazyExtractFrames(this, ipr)
            ip = inputParser;
            addRequired(ip, 'ipr', @isstruct);
            parse(ip, ipr);
            
            fqfps = cell(1, sum(ipr.frames));
            f___ = 0;
            for f = 1:length(ipr.frames)
                if (ipr.frames(f))
                    f___ = f___ + 1;
                    fqfps{f___} = this.frameFileprefix(ipr.dest, f);
                    if (~lexist(this.fourdfpImg(fqfps{f___})))
                        this.buildVisitor.extract_frame_4dfp(ipr.dest, f___, ['-o' fqfps{f___}]);
                    end
                end
            end
            assert(f___ == sum(ipr.frames));
        end
        function [fp,this] = lazyMaskForFrame(this, maskFp, fpm, fpn, m, n)
            if (isempty(maskFp))
                fp = '';
                return
            end
            
            m__ = this.fortranFrameIndices(m);
            n__ = this.fortranFrameIndices(n);
            maskMN = sprintf('%s_%i_%i.4dfp.img', maskFp, m__, n__);
            maskNM = sprintf('%s_%i_%i.4dfp.img', maskFp, n__, m__);
            if (lexist(maskMN, 'file'))
                fp = maskMN;
            else
                if (lexist(maskNM, 'file'))
                    fp = maskNM;
                else
                    [fp,this] = this.maskForFrames(fpm, fpn, maskMN);
                end
            end
            
            %this = this.pushTrash([fp '*.4dfp.*']);
        end
        function [fp,this] = lazyMaskForFrames(this, ipr, fps, m, n)
            ip = inputParser;
            addRequired(ip, 'ipr', @isstruct);
            addRequired(ip, 'fps', @iscell);
            addRequired(ip, 'm',   @isnumeric);
            addRequired(ip, 'n',   @isnumeric);
            parse(ip, ipr, fps, m, n);
            
            m__ = this.fortranFrameIndices(m);
            n__ = this.fortranFrameIndices(n);
            
            if (isempty(ipr.framesMask))
                fp = '';
                return
            end
            maskMN = sprintf('%s_%i_%i.4dfp.img', ipr.framesMask, m__, n__);
            maskNM = sprintf('%s_%i_%i.4dfp.img', ipr.framesMask, n__, m__);
            if (lexist(maskMN, 'file'))
                fp = maskMN;
            else
                if (lexist(maskNM, 'file'))
                    fp = maskNM;
                else
                    [fp,this] = this.maskForFrames(fps{m}, fps{n}, maskMN);
                end
            end            
            
            %this = this.pushTrash([fp '*.4dfp.*']);
        end
        function m         = mprage(this, varargin)
            if (~isempty(this.mprage_))
                m = this.mprage_;
                return
            end
            m = this.sessionData.mprage(varargin{:});
        end
        function this      = msktgenMprage(this, fqfp, atl)
            log = sprintf('msktgenMprage_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
            this.buildVisitor.msktgen2_4dfp(fqfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
        end
        function [fr,this] = nonEmptyFrames(this, varargin)
            %% NONEMPTYFRAMES
            %  @param fqfn is the filename for the native dynamic tracer image;
            %  default is this.nonEmptyFqfilename.
            %  @param named fracThresh is < 1.
            %  @returns fr, a binary vector indicating nonempty frames
            %  @returns this
            
            ip = inputParser;
            addOptional(ip, 'fqfn', this.nonEmptyFqfilename, @(x) lexist(x, 'file'));
            addParameter(ip, 'fracThresh', this.FRAC_FRAME_THRESH, @isnumeric);
            parse(ip, varargin{:});
            
            if (~isempty(this.nonEmptyFrames_))
                fr = this.nonEmptyFrames_;
                return
            end
            
            tr = mlfourd.NumericalNIfTId.load(ip.Results.fqfn);
            tr = tr.volumeSummed;
            maxTr = max(tr.img);
            fr = double(tr.img > ip.Results.fracThresh * maxTr);
            fr = fr';
            this.nonEmptyFrames_ = fr;
        end
        function fqfn      = nonEmptyFqfilename(this)
            %% NONEMPTYFQFILENAME
            %  See also:  mlfourdfp.T4ResolveBuilder0.nonEmptyFrames
            
            sessd = this.sessionData;
            sessd.rnumber = 1;
            if (lexist(sessd.tracerRevision('typ','4dfp.ifh'), 'file'))
                fqfn = sessd.tracerRevision('typ','4dfp.ifh');
            elseif (lexist(sessd.tracerNative('typ','4dfp.ifh'), 'file'))
                fqfn = sessd.tracerNative('typ','4dfp.ifh');
            else 
                error('mlfourdfp:missedRequiredImages', 'T4ResolveBuilder0.nonEmptyFrames');
            end
        end
        function this      = pasteFrames(this, varargin)
            
            ip = inputParser;
            addRequired(ip, 'ipr', @isstruct);
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            ipr = ip.Results.ipr;
            tag = mybasename(ip.Results.tag);
            
            assert(isfield(  ipr, 'dest'));
            assert(ischar(   ipr.dest));
            assert(isfield  (ipr, 'frames'));
            assert(isnumeric(ipr.frames));
            
            pasteList = sprintf('%s_%s_paste.lst', ipr.dest, tag);
            if (lexist(pasteList)); delete(pasteList); end
            
            fid = fopen(pasteList, 'w');
            for f = 1:length(ipr.frames)
                if (ipr.frames(f))
                    fqfp = this.frameFileprefix(ipr.dest, f);
                    fprintf(fid, '%s_%s.4dfp.img\n', fqfp, tag);
                end
            end
            fclose(fid);
            this.buildVisitor.paste_4dfp(pasteList, [ipr.dest '_' tag], 'options', '-a ');
        end
        function dest      = pasteImages(this, fqfps, dest)
            assert(ischar(dest));
            assert(length(fqfps) > 1);
            
            pasteList = [dest '_paste.lst'];
            if (lexist(pasteList))
                delete(pasteList); 
            end
            fid = fopen(pasteList, 'w');
            for idx = 1:length(fqfps)
                fprintf(fid, '%s.4dfp.img\n', fqfps{idx});
            end
            fclose(fid);
            
            this.buildVisitor.paste_4dfp(pasteList, dest, 'options', '-a ');
            %deletefiles(fqfps{:});
        end   
        function fqfp      = sumTimes(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn', @mlfourdfp.FourdfpVisitor.lexist_4dfp);
            addParameter(ip, 'frames', this.frames, @isnumeric);
            parse(ip, varargin{:});
            
            if (sum(ip.Results.frames) < this.readLength(ip.Results.fqfpDyn))
                fqfp = this.sumTimesImagingContext(varargin{:});
            else
                fqfp = this.sumTimesActmapf4dfp(varargin{:});
            end
        end
        function t4s       = t4sInv(this, t4s)
            for t = 1:length(t4s)
                this.buildVisitor.t4_inv(t4s{t});
                t4s{t} = this.buildVisitor.filenameT4inv(t4s{t});
            end
        end
        function t4s       = t4sMul(this, t4s1, t4s2)
            for t = 1:min(length(t4s1), length(t4s2))
                t4s{t} = this.buildVisitor.filenameT4mul(t4s1{t}, t4s2{t}); %#ok<AGROW>
                this.buildVisitor.t4_mul(t4s1{t}, t4s2{t}, t4s{t});                
            end
        end
        function t4s       = t4sMuls(this, varargin)
            nva = length(varargin);
            t4s = varargin{end};
            for v = nva-1:-1:1
                t4s = this.t4sMul(t4s, varargin{v});
            end
        end
        function             t4sImg(this, t4s, source, dest, ref_fdfp)
            assert(iscell(t4s));
            assert(iscell(dest) && all(size(t4s) == size(dest)));
            source    = this.char2cell(source, size(t4s));
            ref_fdfp = this.char2cell(ref_fdfp, size(t4s));
            for f = 1:length(t4s)
                this.buildVisitor.t4img_4dfp(t4s{f}, source{f}, 'out', dest{f}, 'options', ['-O' ref_fdfp{f}]);
            end
        end
        function             touchCompleted(this)
            mlbash(['touch ' this.completedTouchFile]);
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        blurArg_
        buildVisitor_
        frameRegOrdinalFrameF_ = inf
        frameRegOrdinalFrame0_ = 1
        frames_
        mprage_
        nonEmptyFrames_
        product_
        resolveTag_
        sessionData_
        targetFrame_
        trash_ = {};
        xfm_
    end
    
    methods (Access = protected)
        function fp        = blurredFrame(this, fp, varargin)
            ip = inputParser;
            addRequired( ip, 'fp', @ischar);
            addParameter(ip, 'blur', this.blurArg, @isnumeric);
            parse(ip, fp, varargin{:});
            
            fp = this.fileprefixBlurred(ip.Results.fp);
            if (~this.buildVisitor.lexist_4dfp(fp))
                this.buildVisitor.imgblur_4dfp(ip.Results.fp, ip.Results.blur);
            end
        end       
        function c         = char2cell(~, c, sz)
            if (~iscell(c))
                c = cellfun(@(x) c, cell(sz), 'UniformOutput', false);
            end
        end
        function this      = deleteTrash(this)
            for it = 1:length(this.trash_)
               delete(this.trash_{it}); 
            end
            this.trash_ = {};
        end
        function fqfp      = ensureBlurred(this, fqfp)
            %% ENSUREBLURRED
            %  @param fqfp is a fileprefix or cell-array of fileprefixes
            %  @returns fqfp "
            
            if (iscell(fqfp))
                for f = 1:length(fqfp)
                    fqfp{f} = this.ensureBlurred(fqfp{f});
                end
                return
            end
            if (~lstrfind(fqfp, this.fileprefixBlurred('')))
                this.buildVisitor.imgblur_4dfp(fqfp, this.blurArg);
                fqfp = this.fileprefixBlurred(fqfp);
            end
        end
        function fn        = filenameHdr(~, fp)
            fn = [fp '.4dfp.hdr'];
        end
        function fn        = filenameIfh(~, fp)
            fn = [fp '.4dfp.ifh'];
        end
        function fn        = filenameImg(~, fp)
            fn = [fp '.4dfp.img'];
        end
        function fn        = filenameImgRec(~, fp)
            fn = [fp '.4dfp.img.rec'];
        end
        function fp        = iterationFileprefix(this, fp, rnumber)
            %% ITERATIONFILEPREFIX strips this.resolveTag and r[0-9] from fp; replaces r[0-9].
            %  @param fp, e.g., fdgv1r1_resolved
            %  @param rnumber is numeric
            %  @returns fp, e.g., fdgv1r2
            
            assert(ischar(fp));
            assert(isnumeric(rnumber));
            
            idx = regexp(fp, sprintf('_%s$', this.resolveTag));
            if (~isempty(idx))
                fp = fp(1:idx-1);
            end
            if (~isempty(regexp(fp(end-2:end), 'r[0-9]', 'match')))
                fp = fp(1:end-2);
            end
            fp = sprintf('%sr%i', fp, rnumber);
        end
        function fps       = iterationFileprefixes(this, fps, iter)
            assert(iscell(fps));
            for f = 1:length(fps)
                fps{f} = this.iterationFileprefix(fps{f}, iter);
            end
        end
        function msk       = maskBoundaries(this, fqfp)
            if (~this.doMaskBoundaries)
                msk = 'none';
                return
            end
            this.buildVisitor.nifti_4dfp_ng(fqfp);
            ic  = mlfourd.ImagingContext(fqfp);
            ic  = ic.ones;
            ic.noclobber = false;
            ic.saveas([fqfp '_ones']);
            this.buildVisitor.nifti_4dfp_4([fqfp '_ones']);
            msk = this.zeroSlicesOnBoundaries([fqfp '_ones'], 3);
        end
        function [fp,this] = maskForFrames(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fp1', @(x) lexist([x '.4dfp.img']));
            addRequired(ip, 'fp2', @(x) lexist([x '.4dfp.img']));
            addOptional(ip, 'fp',  ['mask_' varargin{1} '_' varargin{2}], @ischar);
            parse(ip, varargin{:});
            
            import mlfourd.* mlfourdfp.*;
            fqfp1 = [ip.Results.fp1 '_g0.1'];
            fqfp2 = [ip.Results.fp2 '_g0.1'];
            bv = this.buildVisitor;
            if (~bv.lexist_4dfp(fqfp1))
                fqfp1 = bv.gauss_4dfp(ip.Results.fp1, 0.1);
            end
            if (~bv.lexist_4dfp(fqfp2))
                fqfp2 = bv.gauss_4dfp(ip.Results.fp2, 0.1);
            end
            ic1 = ImagingContext([fqfp1 '.4dfp.ifh']);
            ic2 = ImagingContext([fqfp2 '.4dfp.ifh']);
            ic  = ImagingContext(ic1.numericalNiftid + ic2.numericalNiftid);
            ic  = ic.thresh(ic.numericalNiftid.dipmax/8);
            ic  = ic.binarized;
            ic.saveas([ip.Results.fp '.4dfp.ifh']);
            fp  = ic.fqfileprefix;
            
            %this = this.pushTrash([fqfp1 '*.4dfp.*']);
            %this = this.pushTrash([fqfp2 '*.4dfp.*']);
        end
        function fqfp      = maskFrameReg(this, ipr)
            if (~isempty(ipr.mprage) && ~isempty(ipr.atlas))
                fqfp = this.fileprefixMsktgen(ipr.dest);
                if (~lexist(this.fourdfpImg(fqfp)))
                    this.msktgenInitial(ipr);
                end
            else
                fqfp = 'none';
            end
        end
        function this      = pushTrash(this, t)
            this.trash_ = [this.trash_ t];
        end
        function fp        = resolvedFileprefix(this, fp, iter)
            assert(ischar(fp));
            assert(isnumeric(iter));
            
            idx = regexp(fp, ['_' this.resolveTag]);
            if (~isempty(idx))
                fp = fp(1:idx-1);
            end
            if (~isempty(regexp(fp(end-2:end), 'r[0-9]', 'match')))
                fp = fp(1:end-2);
            end
            fp = sprintf('%sr%i_%s', fp, iter, this.resolveTag);
        end        
        function this      = t4imgFrames(this, ipr, tag)
            tag = mybasename(tag);
            for f = 1:length(ipr.frames)
                if (ipr.frames(f))
                    frameFp = sprintf('%s_frame%i', ipr.dest, f);
                    this.buildVisitor.t4img_4dfp( ...
                        sprintf('%s_to_%s_t4', frameFp, tag), ...
                        frameFp, ...
                        'out', sprintf('%s_%s', frameFp, tag), ...
                        'options', ['-O' frameFp]);
                end
            end
        end
        function log       = validateForT4Resolve(this, varargin)

            ip = inputParser;
            addParameter(ip, 'dest', 'fdgv1r1', @ischar);
            addParameter(ip, 'frame0', 2, @isnumeric);
            addParameter(ip, 'frameF', 72, @isnumeric);
            parse(ip, varargin{:});
            
            log = mlpipeline.Logger(fullfile(pwd, 'T4ResolveBuilder0.validateForT4Resolve.log'));
            s = sprintf('T4ResolveBuilder0.validateForT4Resolve is starting..........\n');
                fprintf(s);
            log.add(s);
            s = sprintf('work path -> %s\n', pwd);
                fprintf(s);
            log.add(s);
                    
            fdfpVisitor = mlfourdfp.FourdfpVisitor;
            if (isdir('T4'))                
                for f = ip.Results.frame0:ip.Results.frameF            
                    for g = ip.Results.frame0:ip.Results.frameF
                        if (f ~= g)
                            u4 = fullfile('T4', ...
                                fdfpVisitor.filenameT4( ...
                                sprintf('%s_frame%i', ip.Results.dest, g), ...
                                sprintf('%s_frame%i', ip.Results.dest, f)));

                            if (~lexist(u4, 'file'))
                                s = sprintf('Missing file:  %s\n', u4);
                                    fprintf(s);
                                log.add(s);
                            end
                        end
                    end

                    u4 = fullfile('T4', ...
                        fdfpVisitor.filenameT4( ...
                        sprintf('%s_frame%i', ip.Results.dest, f), ...
                        this.resolveTag));

                    if (~lexist(u4, 'file'))
                        s = sprintf('Missing file:  %s\n', u4);
                            fprintf(s);
                        log.add(s);
                    end
                end
            else
                for f = ip.Results.frame0:ip.Results.frameF            
                    for g = ip.Results.frame0:ip.Results.frameF
                        if (f ~= g)
                            u4 = fdfpVisitor.filenameT4( ...
                                sprintf('%s_frame%i', ip.Results.dest, g), ...
                                sprintf('%s_frame%i', ip.Results.dest, f));

                            if (~lexist(u4, 'file'))
                                s = sprintf('Missing file:  %s\n', u4);
                                    fprintf(s);
                                log.add(s);
                            end
                        end
                    end
                end
                for f = ip.Results.frame0:ip.Results.frameF
                    frameFp = sprintf('%s_frame%i', ip.Results.dest, f);
                    if (~lexist([frameFp '.4dfp.ifh'], 'file'))
                        s = sprintf('Missing 4dfp files:  %s\n', frameFp);
                            fprintf(s);
                        log.add(s);
                    end
                end            
            end
            
            log.save;
            diary off;
        end
        function fqfp      = zeroSlicesOnBoundaries(this, fqfp, n)
            N     = this.readLength(fqfp, 3);
            fqfp1 = sprintf('%s_z%ito%i', fqfp,  1, n);
            fqfp2 = sprintf('%s_z%ito%i', fqfp1, N-n+1, N);
            this.buildVisitor.zero_slice_4dfp(fqfp,  'z', 1, n,     fqfp1);
            this.buildVisitor.zero_slice_4dfp(fqfp1, 'z', N-n+1, N, fqfp2);
            fqfp  = fqfp2;            
        end
    end
    
    %% PRIVATE
    
    methods (Access = private)
        function         alignDynamicToAtlas(this, ipr, sumtG, toMprGT4)
            
            mprToAtlT4 = [ipr.mprage '_to_' ipr.atlas '_t4'];
            if (~lexist(mprToAtlT4, 'file'))
                this.msktgenMprage(ipr.mprage, ipr.atlas);
            end
            dynToAtlT4 = [sumtG '_to_' ipr.atlas '_t4'];
            this.buildVisitor.t4_mul(toMprGT4, mprToAtlT4, dynToAtlT4);
            this.buildVisitor.t4img_4dfp(dynToAtlT4, ipr.dest, ...
                'out', [mybasename(ipr.dest), '_on_' mybasename(ipr.atlas)], ...
                'options', ['-O' fullfile(getenv('REFDIR'), ipr.atlas)]);       
        end
        function [sumtG,toMprGT4,log] = ...
                         alignDynamicToMprage(this, ipr)
            
            mprG = this.gaussMprage(ipr);            
            [sumt,sumtG,toMprGT4,toMprGT40] = this.sumTimesDynamic(ipr, mprG);            
            this.buildVisitor.gauss_4dfp(sumt, this.gaussArg, sumtG); 
            
            petMsk = this.maskBoundaries(sumtG);
            log    = sprintf('msktgenInitial_%s.log', datestr(now, 30));
            toMprGT4 = this.buildVisitor.align_multiSpectral( ...
                'dest', mprG, ...
                'source', sumtG, ...
                'destMask', 'none', ...
                'sourceMask', petMsk, ...
                't40', toMprGT40, ...
                't4', toMprGT4, ...
                'log', log);
        end      
        function ipr   = frameRegPar(this, ipr)
            extractedFps = this.lazyExtractFrames(ipr);
            blurredFps   = this.lazyBlurredFrames(extractedFps, ipr);
            parfor m = this.frameRegOrdinalFrame0_:min(length(blurredFps),this.frameRegOrdinalFrameF_)
                for n = 1:length(blurredFps)
                    if (m ~= n) 
                        try
                            maskFp = this.lazyMaskForFrames(ipr, extractedFps, m, n); %#ok<PFBNS>
                            this.buildVisitor.align_2051( ...
                                'dest',       blurredFps{m}, ...
                                'source',     blurredFps{n}, ...
                                'destMask',   maskFp, ...
                                'sourceMask', maskFp, ...
                                't4',         this.buildVisitor.filenameT4(extractedFps{n}, extractedFps{m}), ...
                                't4img_4dfp', false, ...
                                'log',        this.frameRegLog);
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            copyfile( ...
                                this.buildVisitor.transverse_t4, ...
                                this.buildVisitor.filenameT4(extractedFps{n}, extractedFps{m}), 'f');
                            handwarning(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end     
        function ipr   = frameRegSerial(this, ipr)
            extractedFps = this.lazyExtractFrames(ipr);
            blurredFps   = this.lazyBlurredFrames(extractedFps, ipr);
            for m = this.frameRegOrdinalFrame0_:min(length(blurredFps),this.frameRegOrdinalFrameF_)
                for n = 1:length(blurredFps)
                    if (m ~= n) 
                        try
                            t4 = this.buildVisitor.filenameT4(extractedFps{n}, extractedFps{m});
                            if (~lexist(t4))
                                [maskFp,this] = this.lazyMaskForFrames(ipr, extractedFps, m, n);
                                this.buildVisitor.align_2051( ...
                                    'dest',       blurredFps{m}, ...
                                    'source',     blurredFps{n}, ...
                                    'destMask',   maskFp, ...
                                    'sourceMask', maskFp, ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'log',        this.frameRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            copyfile( ...
                                this.buildVisitor.transverse_t4, ...
                                this.buildVisitor.filenameT4(extractedFps{n}, extractedFps{m}), 'f');
                            handwarning(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end
        function ipr   = frameRegSingle(this, ipr)
            try
                extractedFp1 = sprintf('%s_frame%i', ipr.dest, ipr.frame1st);
                extractedFp2 = sprintf('%s_frame%i', ipr.dest, ipr.frame2nd);
                blurredFp1   = this.buildVisitor.fileprefixBlurred(extractedFp1, this.blurArg);
                blurredFp2   = this.buildVisitor.fileprefixBlurred(extractedFp2, this.blurArg);
                mlbash(sprintf('rm -rf %s', this.buildVisitor.filenameT4(extractedFp1, extractedFp2)));
                [maskFp,this] = this.lazyMaskForFrame( ...
                    ipr.framesMask, extractedFp2, extractedFp1, ipr.frame2nd, ipr.frame1st);
                this.buildVisitor.align_2051( ...
                    'dest',       blurredFp2, ...
                    'source',     blurredFp1, ...
                    'destMask',   maskFp, ...
                    'sourceMask', maskFp, ...
                    't4',         this.buildVisitor.filenameT4(extractedFp1, extractedFp2), ...
                    't4img_4dfp', false);
                % t4_resolve requires an idiomatic naming convention for t4 files,
                % based on the names of frame files
                % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4
                
                this.deleteTrash;
                
            catch ME
                copyfile( ...
                    this.buildVisitor.transverse_t4, ...
                    this.buildVisitor.filenameT4(extractedFp1, extractedFp2), 'f');
                handwarning(ME);
            end
        end 
        function fqfp  = gaussMprage(this, ipr)
            if (isempty(ipr.mprage))
                return
            end
            fqfp = this.fileprefixGaussed(ipr.mprage);
            if (~lexist(this.fourdfpImg(fqfp)))
                this.buildVisitor.gauss_4dfp(ipr.mprage, this.gaussArg, fqfp);
            end
        end 
        function         msktgenInitial(this, ipr)
            %% MSKTGENINITIAL 
            %  @params ipr is ip.Results from this.t4ResolveParser
            %  @returns
            %  - blurred MPRAGE named [ipr.mprage '_gNN']
            %  -         summed dynamic imaging named [ipr.dest '_sumt']
            %  - blurred summed dynamic imaging named [ipr.dest '_sumt_gNN']
            %  - t4 files for [ipr.dest '_sumt_gNN'] to ipr.mprage; 
            %                  ipr.mprage to ipr.atlas; 
            %                 [ipr.dest '_sumt_gNN'] to ipr.atlas;
            %                  ipr.dest to ipr.atlas
            %  - blurred summed dynamic imaging named [ipr.dest '_sumt_gNN_on_' ipr.mprage]
            %  -                dynamic imaging named [ipr.dest '_on_' ipr.atlas]
            %  - mask for       dynamic imaging named [ipr.dest '_sumt_gNN_mskt']          
            
            [sumtG,toMprGT4,log] = this.alignDynamicToMprage(ipr);             
                                   this.alignDynamicToAtlas(ipr, sumtG, toMprGT4);
            
            this.buildVisitor.msktgen2_4dfp( ...
                sumtG, this.msktgenThresh, ...
                'options', ['-T' fullfile(getenv('REFDIR'), ipr.atlas)], 'log', log);
        end
        function fqfp  = sumTimesActmapf4dfp(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn');
            addParameter(ip, 'tag', 'sumt', @ischar);
            parse(ip, varargin{:});
            
            fqfp = sprintf('%s_%s', ip.Results.fqfpDyn, ip.Results.tag);
            if (~mlfourdfp.FourdfpVisitor.lexist_4dfp(fqfp))
                this.buildVisitor.actmapf_4dfp( ...
                    sprintf('%i+', this.readLength(ip.Results.fqfpDyn)), ip.Results.fqfpDyn, 'options', ['-a' ip.Results.tag]);
            end
        end
        function [sumt,sumtG,toMprGT4, toMprGT40] = ...
                         sumTimesDynamic(this, ipr, mprG)
                     %% sumTimesDynamic is a KLUDGE that provides housekeeping services to alignDynamicToAtlas.
                     
            mprG  = mybasename(mprG);
            sumt  = mybasename(this.fileprefixSumt(ipr.dest));
            sumtG = mybasename(this.fileprefixGaussed(sumt));
            toMprGT4 = [sumtG '_to_' mprG '_t4'];
            toMprGT40 = ['init_' sumtG '_to_' mprG '_t4'];
            if (lexist(toMprGT4, 'file'))
                delete(toMprGT4);
            end
            
            this.buildVisitor.t4_inv(this.buildVisitor.sagittal_t4, 'out', toMprGT4);
            this.buildVisitor.actmapf_4dfp( ...
                sprintf('%i+', this.readLength(ipr.source)), ipr.dest, 'options', '-asumt', 'normalized', true);
            copyfile(toMprGT4, toMprGT40);
        end
        function fqfp  = sumTimesImagingContext(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn');
            addParameter(ip, 'frames', this.frames);          
            parse(ip, varargin{:});
            fqfp = sprintf('%s_sumt%i', ip.Results.fqfpDyn, sum(ip.Results.frames));
            if (mlfourdfp.FourdfpVisitor.lexist_4dfp(fqfp))
                return
            end
            
            nii = mlfourd.NIfTId([ip.Results.fqfpDyn '.4dfp.ifh']);
            size = nii.size;
            assert(size(4) == length(this.frames));
            img = zeros(size(1:3));
            for f = 1:length(ip.Results.frames)
                if (ip.Results.frames(f))
                    img = img + squeeze(nii.img(:,:,:,f));
                end
            end
            nii.img = img;
            
            nii.saveas([fqfp '.4dfp.ifh']);
        end
    end
    
    %% DEPRECATED
    
    properties (Hidden)
        NIterations = 1
        framesToSkip = 0
    end
    
    methods (Hidden)
        function ipr   = alignDynamicToAtlas0(this, ipr, sumtG, toMprGT4)
            mprToAtlT4 = [ipr.mprage '_to_' ipr.atlas '_t4'];
            if (~lexist(mprToAtlT4, 'file'))
                this.msktgenMprage0(ipr.mprage, ipr.atlas);
            end
            dynToAtlT4 = [sumtG '_to_' ipr.atlas '_t4'];
            this.buildVisitor.t4_mul(toMprGT4, mprToAtlT4, dynToAtlT4);
            this.buildVisitor.t4img_4dfp(dynToAtlT4, ipr.dest, ...
                'out', [ipr.dest, '_on_' ipr.atlas], ...
                'options', ['-O' fullfile(getenv('REFDIR'), ipr.atlas)]);       
        end
        function [sumtG,toMprGT4,log] = ...
                         alignDynamicToMprage0(this, ipr)
            
            mprG = this.gaussMprage(ipr);
            
            [sumt,sumtG,toMprGT4] = this.sumTimesDynamic(ipr, mprG);
            
            this.buildVisitor.gauss_4dfp(sumt, this.gaussArg, sumtG); 
            
            petMsk = this.maskBoundaries(sumtG);
            log    = sprintf('msktgenInitial_%s.log', datestr(now, 30));
            this.buildVisitor.imgreg_4dfp(mprG, 'none', sumtG, petMsk, toMprGT4, 4099, log);
            this.buildVisitor.imgreg_4dfp(mprG, 'none', sumtG, petMsk, toMprGT4, 4099, log);
            this.buildVisitor.imgreg_4dfp(mprG, 'none', sumtG, petMsk, toMprGT4, 3075, log);
            this.buildVisitor.imgreg_4dfp(mprG, 'none', sumtG, petMsk, toMprGT4, 2051, log);            
            this.buildVisitor.imgreg_4dfp(mprG, 'none', sumtG, petMsk, toMprGT4, 2051, log);
            this.buildVisitor.imgreg_4dfp(mprG, 'none', sumtG, petMsk, toMprGT4, 10243, log); 
            this.buildVisitor.t4img_4dfp(toMprGT4, sumtG, 'out', [sumtG '_on_' mprG], 'options', ['-O' mprG]); 
        end
        function fqfps = blurFrames0(this, ipr, fqfps)
            if (~(ipr.blur > 0)); return; end
            for f = 1:length(fqfps)
                if (~isempty(fqfps{f}))
                    this.buildVisitor.imgblur_4dfp(fqfps{f}, this.blurArg);
                    fqfps{f} = this.fileprefixBlurred(fqfps{f});
                end
            end
        end
        function umap  = buildUmapFromCt0(this, ct, umap)
            %% BUILDUMAPFROMCT follows the logic of Lars Couture's ct_umap_4dfp and custom_umap_4dfp.
            %  @param ct   is the (fully-qualified) fileprefix of the CT.
            %  @param umap is the (fully-qualified) fileprefix of the product umap.
            
            assert(lexist(myfilename(ct, '.4dfp.img'), 'file'));
            ct   = myfileprefix(ct);
            assert(ischar(umap) && ~isempty(umap) && ~strcmp(umap, ct));
            umap = myfileprefix(umap);
            ct_  = sprintf('%s_%s', ct, datestr(now, 30));            
            this.buildVisitor.copy_4dfp(ct, ct_);
            %copyfile([ct '.4dfp.ifh'], [ct_ '.4dfp.ifh'], 'f');
            
            this.buildVisitor.scale_4dfp(ct_,  1.0,        'options', '-b-1024');
            this.buildVisitor.scale_4dfp(ct_, -1.0,        'options',         '-ainv');
            this.buildVisitor.scale_4dfp(ct_,  0.079/1326, 'options', '-b0.093 -ap');
            this.buildVisitor.scale_4dfp(ct_,  0.093/1000, 'options', '-b0.093 -an');
            this.buildVisitor.maskimg_4dfp([ct_ '_p'],  ct_,         [ct_ 'p']);
            this.buildVisitor.maskimg_4dfp([ct_ '_n'], [ct_ '_inv'], [ct_ 'n']);
            this.buildVisitor.imgopr_4dfp('a', umap,   [ct_ 'p'],    [ct_ 'n']);
            delete([ct_ '*']);            
        end
        function fqfp  = buildUmapOnSumt0(this)
            ct   = [this.sessionData.ct('typ', 'fqfp') '_on_' mybasename(this.sessionData.mpr('typ', 'fp'))];
            mpr  = this.sessionData.mpr('typ', 'fqfp');
            sumt = this.fileprefixSumt(this.sessionData.fdgNAC('typ', 'fqfp'));
            t40  = this.buildVisitor.t4_inv(fullfile(getenv('RELEASE'), 'S_t4'), 'out', 'S_inv_t4');            
            %t4cm = this.buildVisitor.align_multiSpectral(ct,   mpr, 't40', t40);
            t4sm = this.buildVisitor.align_multiSpectral(sumt, mpr, 't40', t40);
            t4ms = this.buildVisitor.t4_inv(t4sm);
            %t4   = this.buildVisitor.t4_mul(t4cm, t4ms);
            ct   = this.buildVisitor.t4img_4dfp(t4ms, ct, 'options', ['-O' sumt]);
            
            fqfp = this.buildUmapFromCt0(ct, this.sessionData.umapSynth('typ', 'fqfp'));
        end
        function this  = crop(this, ipr)
            if (0 == ipr.crop || ipr.crop == 1)
                if (~lexist(this.filenameIfh(ipr.dest)))
                    dbbash(sprintf('copy_4dfp %s %s', ipr.source, ipr.dest));
                end
                return
            end
                dbbash(sprintf('cropfrac_4dfp %g %s %s', ipr.crop, ipr.source, ipr.dest));
        end
        function fqfps = extractFrames0(this, ipr)
            frameEnd = this.readLength(ipr.source);
            fqfps = cell(1, frameEnd);
            for f = ipr.frame0:ipr.frameF
                fqfps{f} = sprintf('%s_frame%i', ipr.dest, f);
                if (~lexist(this.filenameIfh(fqfps{f})))
                    this.buildVisitor.extract_frame_4dfp(ipr.dest, f);
                end
            end
        end
        function this  = frameReg0(this, ipr, frameFdfps)
            frame0Fdfp = this.filenameIfh(sprintf('%s_frame%i', ipr.dest, ipr.frame0));
            if (~lexist(frame0Fdfp))
                error('mlfourdfp:missingFile', 'T4ResolveBuilder0.frameReg0 could not find %s', frame0Fdfp);
            end            
            blurredFdfps = this.blurFrames0(ipr, frameFdfps);
            log = mlpipeline.Logger( ...
                loggerFilename(ipr.dest, 'func', 'frameReg0'), this);            
            mask = this.frameRegMask(ipr);
            for m = ipr.frame0:length(blurredFdfps)
                for n = ipr.frame0:length(blurredFdfps)
                    if (m ~= n)
                        t4file = sprintf('%s_to_%s_t4', frameFdfps{n}, frameFdfps{m});
                        this.buildVisitor.imgreg_4dfp(blurredFdfps{m}, mask, blurredFdfps{n}, mask, t4file, 2051, log.fqfilename);
                        this.buildVisitor.imgreg_4dfp(blurredFdfps{m}, mask, blurredFdfps{n}, mask, t4file, 2051, log.fqfilename);
                    end
                end
            end
        end
        function fqfp  = frameRegMask(this, ipr) %#ok<INUSD>
            fqfp = 'none';
            
            %fqfp = this.fileprefixMskt(this.fileprefixGaussed(this.fileprefixSumt(ipr.dest)));
            %if (~lexist(this.filenameIfh(fqfp)))
            %    error('mlfourdfp:missingFile', 'T4ResolveBuilder0.frameReg0 could not find %s', maskfn);
            %end
        end
        function this  = msktgenInitial0(this, ipr)
            %% MSKTGENINITIAL 
            %  @params ipr is ip.Results from this.t4ResolveParser
            %  @returns
            %  - blurred MPRAGE named [ipr.mprage '_gNN']
            %  -         summed dynamic imaging named [ipr.dest '_sumt']
            %  - blurred summed dynamic imaging named [ipr.dest '_sumt_gNN']
            %  - t4 files for [ipr.dest '_sumt_gNN'] to ipr.mprage; 
            %                  ipr.mprage to ipr.atlas; 
            %                 [ipr.dest '_sumt_gNN'] to ipr.atlas;
            %                  ipr.dest to ipr.atlas
            %  - blurred summed dynamic imaging named [ipr.dest '_sumt_gNN_on_' ipr.mprage]
            %  -                dynamic imaging named [ipr.dest '_on_' ipr.atlas]
            %  - mask for       dynamic imaging named [ipr.dest '_sumt_gNN_mskt']          
            
            [sumtG,toMprGT4,log] = this.alignDynamicToMprage0(ipr);             
            this                 = this.alignDynamicToAtlas0(ipr, sumtG, toMprGT4);    
            
            this.buildVisitor.msktgen2_4dfp(sumtG, 20, 'options', ['-T' fullfile(getenv('REFDIR'), ipr.atlas)], 'log', log);
        end
        function this  = msktgenMprage0(this, fqfp, atl)
            log = sprintf('msktgenMprage_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
            this.buildVisitor.msktgen2_4dfp(fqfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
        end
        function this  = pasteFrames0(this, varargin)
            %% PASTEFRAMES uses a tag for frame-files
            
            ip = inputParser;
            addRequired(ip, 'ipr', @isstruct);
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            assert(isfield(  ip.Results.ipr, 'dest'));
            assert(ischar(   ip.Results.ipr.dest));
            assert(isfield  (ip.Results.ipr, 'frames'));
            assert(isnumeric(ip.Results.ipr.frames));
            
            if (isempty(ip.Results.tag))
                this = this.pasteFrames__(ip.Results.ipr);
                return
            end
            this = this.pasteFramesTagged__(ip.Results.ipr, ip.Results.tag);
        end
        function this  = pasteFrames__(this, ipr)
            pasteList = sprintf('%s_paste.lst', ipr.dest);
            if (lexist(pasteList)); delete(pasteList); end

            fid = fopen(pasteList, 'w');
            for f = 1:length(ipr.frames)
                if (~ipr.frames(f))
                    ipr0 = ipr; 
                    ipr0.rnumber = max(ipr.rnumber-1, 1);
                    fqfp = this.lazyExtractFrame(ipr0.dest, f);
                    this.buildVisitor.move_4dfp(fqfp, [fqfp '_' this.resolveTag]);
                else
                    fqfp = this.frameFileprefix(ipr.dest, f);
                end
                fprintf(fid, '%s.4dfp.img\n', fqfp);
            end
            fclose(fid);            
            this.buildVisitor.paste_4dfp(pasteList, ipr.dest, 'options', '-a ');
        end
        function this  = pasteFramesTagged__(this, ipr, tag)
            tag = mybasename(tag);
            pasteList = sprintf('%s_%s_paste.lst', ipr.dest, tag);
            if (lexist(pasteList)); delete(pasteList); end
            
            fid = fopen(pasteList, 'w');
            for f = 1:length(ipr.frames)
                if (~ipr.frames(f))
                    ipr0 = ipr; 
                    ipr0.rnumber = max(ipr.rnumber-1, 1);
                    fqfp = this.lazyExtractFrame(ipr0.dest, f);
                    this.buildVisitor.move_4dfp(fqfp, [fqfp '_' this.resolveTag]);
                else
                    fqfp = this.frameFileprefix(ipr.dest, f);
                end
                fprintf(fid, '%s_%s.4dfp.img\n', fqfp, tag);
            end
            fclose(fid);
            this.buildVisitor.paste_4dfp(pasteList, [ipr.dest '_' tag], 'options', '-a ');
        end
        function this  = t4ResolveAndPaste0(this, ipr)
            dest = mybasename(ipr.dest);
            imgFns = '';
            for f = 1:length(ipr.frames)
                imgFns = [imgFns this.filenameImg(sprintf(' %s_frame%i', dest, f))]; %#ok<AGROW>
            end
            
            log = mlpipeline.Logger( ...
                loggerFilename(dest, 'func', 't4ResolveAndPaste'), this);
            
            resolveSuff = 'resolved';
            fprintf('t4_resolve -v -m -s -o%s %s &> %s', resolveSuff, imgFns, log.fqfilename);
            this.buildVisitor.t4_resolve(resolveSuff, imgFns, 'options', '-v -m -s', 'log', log.fqfilename);
            this.t4imgFrames(ipr, resolveSuff);
            this.pasteFrames(ipr, resolveSuff);

            resolved = sprintf('%s_%s_%s', resolveSuff, dest, datestr(now, 30)); 
            copyfile('resolved.mat0', [resolved '.mat0']);
            copyfile('resolved.sub',  [resolved '.sub']);
        end
        function this  = t4ResolveIterative(this, source, dest, mpr)
            source  = mybasename(source);
            dest  = mybasename(dest);
            frame0 = this.framesToSkip+1; %#ok<PROPLC>
            frameF = this.readLength(source); %#ok<PROPLC>
            
            iterations = 1;            
            this = this.t4ResolveIterate( ...
                source, ...
                sprintf('%sr%i', dest, iterations), ...
                'mprage', mpr, ...
                'frame0', frame0, ...
                'crop', this.firstCrop); %#ok<PROPLC>
            
            while (iterations < this.NIterations)
                iterations = iterations + 1;
                if (iterations < 3)
                    this = this.t4ResolveIterate( ...
                        sprintf('%sr%i_%s', dest, iterations-1, frame0, frameF, this.resolveTag), ...
                        sprintf('%sr%i', dest, iterations), ...
                        'mprage', mpr, ...
                        'frame0', this.framesToSkip+1, ...
                        'frameF', frameF-this.framesToSkip); %#ok<PROPLC>
                else
                    this = this.t4ResolveIterate( ...
                        sprintf('%sr%i_%s', dest, iterations-1, 1, frameF-this.framesToSkip, this.resolveTag), ...
                        sprintf('%sr%i', dest, iterations), ...
                        'mprage', mpr, ...
                        'frame0', this.framesToSkip+1, ...
                        'frameF', frameF-this.framesToSkip); %#ok<PROPLC>
                end
            end
            this.product_ = [this.product_ this.filenameImg(sprintf('%sr%i', dest, iterations))];
        end
        function this  = t4ResolveIterate(this, varargin)
            ip = this.t4ResolveParser(varargin{:});          
            dest = ip.Results.dest;
            
            if (~lexist(this.filenameIfh(dest)))
                this = this.crop(ip.Results);
            end
            if (~strcmp('none', this.frameRegMask(ip.Results)) && ...
                ~lexist(this.filenameIfh(this.fileprefixMskt(this.fileprefixGaussed(this.fileprefixSumt(dest)))), 'file'))
                this = this.msktgenInitial0(ip.Results);
            end
            frameFps = this.extractFrames0(ip.Results);
            this = this.frameReg0(ip.Results, frameFps);
            this = this.t4ResolveAndPaste0(ip.Results); 
            this = this.teardownT4ResolveIteration0(ip.Results);
            this.buildVisitor.imgblur_4dfp(dest, this.blurArg);
        end  
        function ip    = t4ResolveParser(this, varargin)
            ip = inputParser;
            addRequired( ip, 'source',  @(x) lexist(this.filenameIfh(x), 'file'));
            addRequired( ip, 'dest',    @ischar);
            addParameter(ip, 'mprage',  this.sessionData.mpr('typ', 'fqfn'), ...
                                        @(x) lexist(this.filenameIfh(x), 'file'));
            addParameter(ip, 'frame0',  1, ...
                                        @isnumeric);
            addParameter(ip, 'frameF',  this.readLength(varargin{1}), ...
                                        @isnumeric);
            addParameter(ip, 'crop',    1, ...
                                        @isnumeric);
            addParameter(ip, 'atlas',   'TRIO_Y_NDC', ...
                                        @(x) lexist(fullfile(getenv('REFDIR'), this.filenameIfh(x))));
            addParameter(ip, 'blur',    this.blurArg, ...
                                        @isnumeric);
            addParameter(ip, 'rnumber', ...
                                        0, @isnumeric );
            parse(ip, varargin{:});  
        end   
        function this  = teardownT4ResolveIteration0(this, ipr)
            if (this.keepForensics); return; end
            
            name = ipr.dest;
            for f = 1:length(ipr.frames)
                delete(sprintf('%s_frame%i.4dfp.*', name, f));
            end
            for f = 1:length(ipr.frames)
                delete(sprintf('%s_frame%i_b*.4dfp.*', name, f));
                delete(sprintf('%s_frame%i_.4dfp.*', name, f, this.resolveTag));
            end
            if (~isdir('T4'))
                mkdir('T4');
                movefiles('*_t4', 'T4');
            end
            this.buildVisitor.imgblur_4dfp(ipr.dest, this.blurArg);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

 