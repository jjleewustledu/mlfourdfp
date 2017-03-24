classdef (Abstract) AbstractT4ResolveBuilder < mlfourdfp.IT4ResolveBuilder
	%% ABSTRACTT4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Jan-2017 22:33:57
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties
        doMaskBoundaries = false
        finished
        keepForensics
        logger % TODO
        mprToAtlT4
        msktgenThresh = 0
        NRevisions 
        
        imageRegLog
        resolveLog
    end
    
    properties (Dependent)
        blurArg
        buildVisitor        
        gaussArg
        imageComposite
        indicesLogical
        indexOfReference
        petBlur
 		product
        referenceImage
        referenceWeight
        resolveTag
        sessionData
        sourceImage
        sourceWeight
        studyData
        theImages
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
        function v    = get.buildVisitor(this)
            v = this.buildVisitor_;
        end
        function this = set.buildVisitor(this, v)
            assert(isa(v, 'mlfourdfp.FourdfpVisitor'));
            this.buildVisitor_ = v;
        end
        function g    = get.gaussArg(this)
            %g = 1.1;
            g = this.F_HALF_x_FWHM/this.blurArg;
        end
        function g    = get.imageComposite(this)
            g = this.imageComposite_;
        end
        function this = set.imageComposite(this, s)
            assert(isa(s, 'mlfourdfp.IImageComposite'));
            this.imageComposite_ = s;
        end
        function g    = get.indexOfReference(this)
            g = this.imageComposite.indexOfReference;
        end
        function this = set.indexOfReference(this, s)
            this.imageComposite_.indexOfReference = s;
        end
        function g    = get.indicesLogical(this)
            g = this.imageComposite.indicesLogical;
        end
        function this = set.indicesLogical(this, s)
            this.imageComposite.indicesLogical = s;
        end
        function g    = get.petBlur(this)
            g = this.sessionData.petBlur;
        end
        function prod = get.product(this)
            prod = this.product_;
        end
        function rw   = get.referenceWeight(~)
            rw = [];
        end
        function ri   = get.referenceImage(this)
            ri = this.imageComposite.referenceImage;
        end
        function g    = get.resolveTag(this)
            g = this.sessionData.resolveTag;
        end
        function this = set.resolveTag(this, s)
            assert(ischar(s));
            this.sessionData_.resolveTag = s;
        end
        function g    = get.sessionData(this)
            g = this.sessionData_;
        end
        function this = set.sessionData(this, s)
            assert(isa(s, 'mlpipeline.SessionData'));
            this.sessionData_ = s;
        end
        function si   = get.sourceImage(this)
            si = this.imageComposite.sourceImage;
        end
        function sw   = get.sourceWeight(~)
            sw = [];
        end
        function g    = get.studyData(this)
            g = this.sessionData.studyData;
        end
        function g    = get.theImages(this)
            g = this.imageComposite.theImages;
        end
        function this = set.theImages(this, s)
            this.imageComposite_.theImages = s;
        end
    end
    
    methods (Static)
        function      diaryv(name)
            assert(ischar(name));
            if (~isempty(getenv('PRINTV')))
                diary(sprintf('mlfourdfp_AbstractT4ResolveBuilder%s.txt', name)); 
            end
        end
        function fn = fourdfpHdr(fp)
            fn = [fp '.4dfp.hdr'];
        end
        function fn = fourdfpIfh(fp)
            fn = [fp '.4dfp.ifh'];
        end
        function fn = fourdfpImg(fp)
            fn = [fp '.4dfp.img'];
        end
        function fn = fourdfpImgRec(fp)
            fn = [fp '.4dfp.img.rec'];
        end
        function f  = frameNumber(str, offset)
            names = regexp(str, '\w+(-|_)(F|f)rame(?<f>\d+)', 'names');
            f = str2double(names.f) + offset;
        end
        function      printv(varargin)
            if (~isempty(getenv('PRINTV')))
                fprintf('mlfourdfp.AbstractT4ResolveBuilder.');
                fprintf(varargin{:});
                fprintf('\n');
            end
        end
    end
    
	methods 
		function this = AbstractT4ResolveBuilder(varargin)
 			%% ABSTRACTT4RESOLVEBUILDER
 			%  @param named blurArg
 			%  @param named buildVisitor
            %  @param named sessionData
            %  @param named indicesLogical
            %  @param named indexOfReference
            %  @param named NRevisions
            %  @param named keepForensics
            %  @param named resolveTag
 			
            if (0 == nargin); return; end
            
            %% invoke copy-ctor
            
            if (1 == nargin && isa(varargin{1}, 'mlfourdfp.AbstractT4ResolveBuilder'))
                aCopy = varargin{1};
                
                %% properties
                this.doMaskBoundaries = aCopy.doMaskBoundaries;
                this.finished = aCopy.finished;
                this.keepForensics = aCopy.keepForensics;
                this.logger = aCopy.logger;
                this.mprToAtlT4 = aCopy.mprToAtlT4;
                this.msktgenThresh = aCopy.msktgenThresh;
                this.NRevisions = aCopy.NRevisions;
                
                this.imageRegLog = aCopy.imageRegLog;
                this.resolveLog = aCopy.resolveLog;
                
                %% properties (Access = protected)
                this.blurArg_ = aCopy.blurArg_;
                this.buildVisitor_ = aCopy.buildVisitor_;
                this.imageComposite_ = aCopy.imageComposite_;
                this.product_ = aCopy.product_;
                this.sessionData_ = aCopy.sessionData_;
                this.trash_ = aCopy.trash_;
                this.xfm_ = aCopy.xfm_; 
                return
            end
            
            %% manage parameters 
            
            import mlfourdfp.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'buildVisitor',     FourdfpVisitor, @(x) isa(x, 'mlfourdfp.FourdfpVisitor'));
            addParameter(ip, 'sessionData',      [],             @(x) isa(x, 'mlpipeline.ISessionData'));
            addParameter(ip, 'NRevisions',       2,              @isnumeric);
            addParameter(ip, 'keepForensics',    true,           @islogical);
            %addParameter(ip, 'resolveTag',       '',             @ischar);
            addParameter(ip, 'theImages',        {},             @(x) iscell(x) || ischar(x));
            parse(ip, varargin{:});
            
            this.buildVisitor_ = ip.Results.buildVisitor;
            this.sessionData_ = ip.Results.sessionData;
            this.sessionData_.builder = this;  
            this.NRevisions = ip.Results.NRevisions;
            this.keepForensics = ip.Results.keepForensics;
            %if (~isempty(ip.Results.resolveTag))
            %    this.resolveTag = ip.Results.resolveTag;
            %end
            this.theImages = FourdfpVisitor.ensureSafeOn(ip.Results.theImages);
            
            this = this.mpr2atl;
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
        function        teardownRevision(this)
            if (this.keepForensics); return; end
            
            this.teardownLogs;
            this.teardownT4s;
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
        function tf   = isfinished(this, varargin)
            tf = this.finished.isfinished;
        end
        
        %% UTILITY
        
        function a     = atlas(this, varargin) 
            a = this.sessionData.atlas(varargin{:});
        end
        function [s,r] = ensureNifti(this, varargin)
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
                  'T4ResolveBuilder.ensureNifti could not find files of form %s', ip.Results.filename);            
        end
        function [s,r] = ensure4dfp(this, varargin)
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
                  'T4ResolveBuilder.ensureNifti could not find files of form %s', ip.Results.filename);     
        end
        function ipr   = expandBlurs(this, ipr)
            assert(isstruct(ipr));
            assert(isfield( ipr, 'sourceBlur'));
            assert(isfield( ipr, 'destBlur'));
            ipr.sourceBlur = ipr.sourceBlur .* this.indicesLogical;
            ipr.destBlur   = ipr.destBlur   .* this.indicesLogical;
        end
        function fn    = filenameHdr(~, fp)
            fn = [fp '.4dfp.hdr'];
        end
        function fn    = filenameIfh(~, fp)
            fn = [fp '.4dfp.ifh'];
        end
        function fn    = filenameImg(~, fp)
            fn = [fp '.4dfp.img'];
        end
        function fn    = filenameImgRec(~, fp)
            fn = [fp '.4dfp.img.rec'];
        end  
        function fqfp  = fileprefixBlurred(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @ischar);
            addOptional(ip, 'blur', this.blurArg, @isnumeric);
            parse(ip, varargin{:});
            fqfp = this.buildVisitor.fileprefixBlurred(ip.Results.fqfp, ip.Results.blur);
        end        
        function fqfp  = fileprefixGaussed(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @ischar);
            addOptional(ip, 'gauss', this.gaussArg, @isnumeric);
            parse(ip, varargin{:});
            fqfp = this.buildVisitor.fileprefixGaussed(ip.Results.fqfp, ip.Results.gauss);
        end
        function fqfp  = fileprefixMsk(~, fqfp)
            fqfp = [fqfp '_msk'];
        end
        function fqfp  = fileprefixMskt(~, fqfp)
            fqfp = [fqfp '_mskt'];
        end
        function fqfp  = fileprefixMsktgen(this, fqfp)
            fqfp = this.fileprefixMskt( ...
                this.fileprefixGaussed(this.fileprefixSumt(fqfp)));
        end
        function fp    = fileprefixResolved(this, fp, rnumber)
            fp = sprintf('%s_%s', this.fileprefixRevision(fp, rnumber), this.resolveTag);
        end   
        function fps   = fileprefixesResolved(this, fps, rnumber)
            %% FILEPREFIXREVISION strips and replaces r[0-9]; adds this.resolveTag.
            %  @param fp, e.g., fdgv1r1_resolved
            %  @param rnumber is numeric
            %  @returns fp, e.g., fdgv1r2_resolved
            %  See also:  mlfourdfp.T4ResolveBuilder.fileprefixRevision
            
            assert(iscell(fps));
            for f = 1:length(fps)
                fps{f} = this.fileprefixResolved(fps{f}, rnumber);
            end
        end 
        function fp    = fileprefixRevision(this, fp, rnumber)
            %% FILEPREFIXREVISION strips this.resolveTag and r[0-9] from fp; replaces r[0-9].
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
        function fps   = fileprefixesRevision(this, fps, rnumber)
            assert(iscell(fps));
            for f = 1:length(fps)
                fps{f} = this.fileprefixRevision(fps{f}, rnumber);
            end
        end
        function fqfp  = fileprefixSumt(this, fqfp)
            assert(ischar(fqfp));
            fqfp = sprintf('%s_sumt%i', fqfp, sum(this.indicesLogical));
        end
        function fqfp  = fileprefixTmp(this, varargin)
            ip = inputParser;
            addParameter(ip, 'path', this.sessionData.sessionPath, @isdir);
            parse(ip, varargin{:});
            
            fqfp = fullfile(ip.Results.path, sprintf('tmp_%s', datestr(now, 30)));            
        end    
        function this  = img2atl(this, fqfp)
            sd = this.sessionData;
            imgToAtlT4 = [fqfp '_to_' sd.atlas('typ', 'fp') '_t4'];
            if (~lexist(imgToAtlT4, 'file'))
                pwd_ = pushd(fileparts(fqfp));
                this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-crossmodal -T' sd.atlas('typ', 'fqfp')]);
                popd(pwd_);
            end
            if (~lexist(mybasename(imgToAtlT4), 'file'))
                this.buildVisitor.lns(imgToAtlT4);
            end
        end
        function fqfps = lazyBlurImages(this, ipr)
            %% LAZYBLURIMAGES uses specifiers in ipr; will not replace any existing indicesLogical
            %  @param ipr is a struct
            
            fqfps = {};
            for il = 1:length(this.indicesLogical)
                if (this.indicesLogical(il))
                    ipr.currentIndex = il;
                    fqfps = [fqfps this.lazyBlurImage(ipr, ipr.destBlur(il))]; %#ok<AGROW>
                end
            end
            assert(length(fqfps) == sum(this.indicesLogical));
        end    
        function msk   = maskBoundaries(this, fqfp)
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
        function fp    = maskForImages(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fp1', @(x) lexist([x '.4dfp.ifh']));
            addRequired(ip, 'fp2', @(x) lexist([x '.4dfp.ifh']));
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
        end        
        function this  = mpr2atl(this)
            sd = this.sessionData;
            this.mprToAtlT4 = [sd.mpr('typ', 'fqfp') '_to_' sd.atlas('typ', 'fp') '_t4'];
            if (~lexist(this.mprToAtlT4, 'file'))
                pwd_ = pushd(sd.mpr('typ', 'path'));
                this.buildVisitor.mpr2atl_4dfp(sd.mpr('typ', 'fp'), 'options', ['-T' sd.atlas('typ', 'fqfp')]);
                popd(pwd_);
            end
            if (~lexist(mybasename(this.mprToAtlT4), 'file'))
                this.buildVisitor.lns(this.mprToAtlT4);
            end
        end
        function this  = msktgenImg(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @lexist_4dfp);
            addOptional(ip, 'atl', fullfile(getenv('REFDIR'), 'TRIO_Y_NDC'), @lexist_4dfp);
            parse(ip, varargin{:});
            fqfp = ip.Results.fqfp;
            atl  = ip.Results.atl;
            
            log = sprintf('msktgenImg_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-crossmodal -T' atl], 'log', log);
            this.buildVisitor.msktgen_4dfp(fqfp, 'options', ['-T' atl], 'log', log);
        end 
        function this  = msktgenMprage(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @lexist_4dfp);
            addOptional(ip, 'atl', fullfile(getenv('REFDIR'), 'TRIO_Y_NDC'), @lexist_4dfp);
            parse(ip, varargin{:});
            fqfp = ip.Results.fqfp;
            atl  = ip.Results.atl;
            
            log = sprintf('msktgenMprage_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fqfp, 'options', ['-T' atl], 'log', log);
            this.buildVisitor.msktgen_4dfp(fqfp, 'options', ['-T' atl], 'log', log);
        end        
        function fqfp  = sumTimes(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn', @lexist_4dfp);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @isnumeric);
            parse(ip, varargin{:});
            
            if (sum(ip.Results.indicesLogical) < this.imageComposite.length)
                fqfp = this.sumTimesImagingContext(varargin{:});
            else
                fqfp = this.sumTimesActmapf4dfp(varargin{:});
            end
        end
        function out   = t4img_4dfp(this, varargin)
            switch (this.NRevisions)
                case 1
                    out = this.t4img_4dfpr1(varargin{:});
                case 2
                    out = this.t4img_4dfpr2(varargin{:});
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', ...
                          'CompositeT4ResolveBuilder.t4img_4dfp.NRevisions->%i', this.NRevisions);
            end
        end
        function out   = t4img_4dfpr1(this, varargin)
            ip = inputParser;
            addRequired(ip, 'src0', @(x) lstrfind(x, this.theImages));
            addRequired(ip, 'src',  @lexist_4dfp);
            addOptional(ip, 'ref', this.theImages{1}, @lexist_4dfp);
            addParameter(ip, 'opts', '', @ischar);
            parse(ip, varargin{:});            
            src0 = ip.Results.src0;
            src  = ip.Results.src;
            ref  = ip.Results.ref;
            opts = ip.Results.opts;
            
            t4  = sprintf('%sr1_to_%s_t4', src0, this.resolveTag);
            in  = src;
            out = sprintf('%sr1_%s', this.clipLastRevisionMarking(src), this.resolveTag);
            this.buildVisitor.t4img_4dfp(t4, in, 'out', out, 'options', ['-O' ref 'r1_' this.resolveTag ' ' opts]);
        end
        function out   = t4img_4dfpr2(this, varargin)
            ip = inputParser;
            addRequired(ip, 'src0', @(x) lstrfind(x, this.theImages));
            addRequired(ip, 'src',  @lexist_4dfp);
            addOptional(ip, 'ref', this.theImages{1}, @lexist_4dfp);
            addParameter(ip, 'opts', '', @ischar);
            parse(ip, varargin{:});            
            src0 = ip.Results.src0;
            src  = ip.Results.src;
            ref  = ip.Results.ref;
            opts = ip.Results.opts;
            
            t4  = sprintf('%sr1_to_%s_t4', src0, this.resolveTag);
            in  = src;
            out = sprintf('%sr1_%s', this.clipLastRevisionMarking(src), this.resolveTag);
            this.buildVisitor.t4img_4dfp(t4, in, 'out', out, 'options', ['-O' ref 'r1_' this.resolveTag ' ' opts]);
            
            t4  = sprintf('%sr2_to_%s_t4', src0, this.resolveTag);
            in  = sprintf('%sr1_%s', this.clipLastRevisionMarking(src), this.resolveTag);
            out = sprintf('%sr2_%s', this.clipLastRevisionMarking(src), this.resolveTag);
            this.buildVisitor.t4img_4dfp(t4, in, 'out', out, 'options', ['-O' ref 'r2_' this.resolveTag ' ' opts]);
        end
 	end 

    %% PROTECTED
    
    properties (Access = protected)
        blurArg_
        buildVisitor_        
        imageComposite_
        product_
        sessionData_
        trash_ = {};
        xfm_
    end
    
    methods (Access = protected)
        function fp   = clipLastRevisionMarking(~, fp)
            pos = regexp(fp, 'r\d$', 'ONCE');
            if (~isempty(pos))
                fp = fp(1:pos-1);
            end
        end
        function this = deleteTrash(this)
            for it = 1:length(this.trash_)
               delete(this.trash_{it}); 
            end
            this.trash_ = {};
        end
        function this = pushTrash(this, t)
            this.trash_ = [this.trash_ t];
        end
        function fqfp = sumTimesActmapf4dfp(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn');
            addParameter(ip, 'tag', 'sumt', @ischar);
            parse(ip, varargin{:});
            
            fqfp = sprintf('%s_%s', ip.Results.fqfpDyn, ip.Results.tag);
            if (~lexist_4dfp(fqfp))
                this.buildVisitor.actmapf_4dfp( ...
                    sprintf('%i+', this.imageComposite.length, ip.Results.fqfpDyn, 'options', ['-a' ip.Results.tag]));
            end
        end
        function fqfp = sumTimesImagingContext(this, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfpDyn');
            addParameter(ip, 'indicesLogical', this.indicesLogical);          
            parse(ip, varargin{:});
            fqfp = sprintf('%s_sumt%i', ip.Results.fqfpDyn, sum(ip.Results.indicesLogical));
            if (lexist_4dfp(fqfp))
                return
            end
            
            nii = mlfourd.NIfTId([ip.Results.fqfpDyn '.4dfp.ifh']);
            size = nii.size;
            assert(size(4) == length(this.indicesLogical));
            img = zeros(size(1:3));
            for f = 1:length(ip.Results.indicesLogical)
                if (ip.Results.indicesLogical(f))
                    img = img + squeeze(nii.img(:,:,:,f));
                end
            end
            nii.img = img;
            
            nii.saveas([fqfp '.4dfp.ifh']);
        end
        function fqfp = zeroSlicesOnBoundaries(this, fqfp, n)
            N     = this.imageComposite.length;
            fqfp1 = sprintf('%s_z%ito%i', fqfp,  1, n);
            fqfp2 = sprintf('%s_z%ito%i', fqfp1, N-n+1, N);
            this.buildVisitor.zero_slice_4dfp(fqfp,  'z', 1, n,     fqfp1);
            this.buildVisitor.zero_slice_4dfp(fqfp1, 'z', N-n+1, N, fqfp2);
            fqfp  = fqfp2;            
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

