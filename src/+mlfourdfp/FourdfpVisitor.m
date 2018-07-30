classdef FourdfpVisitor 
	%% FOURDFPVISITOR  

	%  $Revision$
 	%  was created 01-Mar-2016 18:43:03
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.  Copyright 2017-2018 John Joowon Lee.
 	

	properties (Constant)
        ASSERT_PLATFORM = false
 		FOURDFP_HOSTS = ...
            {'william' 'acom' 'mahler' 'maulinux1' 'pascal' 'linux5' 'bruckner' 'wagner' 'cluster'}
        UNSAFE_FILEPREFIXES = {'_on_' '+' '.a2009s'} % '_op_' removed 2018feb8
        SUPPORTED_EXT = {'.4dfp.ifh' '.4dfp.hdr' '.4dfp.img' '.4dfp.img.rec'}
    end
    
    properties 
        coronal_t4      = fullfile(getenv('RELEASE'), 'C_t4')
        sagittal_t4     = fullfile(getenv('RELEASE'), 'S_t4')
        sagittal_inv_t4 = fullfile(getenv('HOME'), 'Local', 'bin', 'Sinv_t4')
        transverse_t4   = fullfile(getenv('RELEASE'), 'T_t4')
    end
    
    methods (Static)
        function fqfp  = backupn(varargin)
            %% BACKUPN makes backup copies of 4dfp files, keeping as many as n backup versions.
            %  @param fqfp is a 4dfp fileprefix
            %  @param n is an integer
            %  @throws error mlfourdfp:backupRejected
            
            import mlfourdfp.*;
            ip = inputParser;
            addRequired( ip, 'fqfp', @FourdfpVisitor.lexist_4dfp);
            addOptional( ip, 'n', 1, @isnumeric);
            addParameter(ip, 'tag', '_backup', @ischar);
            parse(ip, varargin{:});
            
            dt = mlsystem.DirTool([ip.Results.fqfp ip.Results.tag '*.4dfp.img']);
            if (length(dt.fqfns) < ip.Results.n)
                this = FourdfpVisitor;
                fqfp = sprintf('%s%sD%s', ip.Results.fqfp, ip.Results.tag, datestr(now, 30));
                this.copy_4dfp(ip.Results.fqfp, fqfp);
                return
            end
            error('mlfourdfp:backupRejected', 'FourdfpVisitor.backupn.ip.Results.n -> %i', ip.Results.n);
        end
        function         copyfile_4dfp(varargin)
            v = varargin;
            xs = mlfourdfp.FourdfpVisitor.SUPPORTED_EXT;
            v = cellfun(@(x) myfileprefix(x), v, 'UniformOutput', false);
            if (nargin > 1 && isdir(v{2}))
                for ix = 1:length(xs)
                    copyfile([v{1} xs{ix}], v{2})
                end
                return
            end
            for ix = 1:length(xs)
                fns = cellfun(@(x) [x xs{ix}], v, 'UniformOutput', false);
                copyfile(fns{:});
            end
        end
        function         copyfilef_4dfp(varargin)
            v = varargin;
            xs = mlfourdfp.FourdfpVisitor.SUPPORTED_EXT;
            v = cellfun(@(x) myfileprefix(x), v, 'UniformOutput', false);
            if (nargin > 1 && isdir(v{2}))
                for ix = 1:length(xs)
                    copyfile([v{1} xs{ix}], v{2}, 'f')
                end
                return
            end
            for ix = 1:length(xs)
                fns = cellfun(@(x) [x xs{ix}], v, 'UniformOutput', false);
                copyfile(fns{:}, 'f');
            end
        end
        function         movefile_4dfp(varargin)
            v = varargin;
            xs = mlfourdfp.FourdfpVisitor.SUPPORTED_EXT;
            v = cellfun(@(x) myfileprefix(x), v, 'UniformOutput', false);
            if (nargin > 1 && isdir(v{2}))
                for ix = 1:length(xs)
                    movefile([v{1} xs{ix}], v{2})
                end
                return
            end
            for ix = 1:length(xs)
                fns = cellfun(@(x) [x xs{ix}], v, 'UniformOutput', false);
                movefile(fns{:});
            end
        end
        function         movefilef_4dfp(varargin)
            v = varargin;
            xs = mlfourdfp.FourdfpVisitor.SUPPORTED_EXT;
            v = cellfun(@(x) myfileprefix(x), v, 'UniformOutput', false);
            if (nargin > 1 && isdir(v{2}))
                for ix = 1:length(xs)
                    movefile([v{1} xs{ix}], v{2}, 'f')
                end
                return
            end
            for ix = 1:length(xs)
                fns = cellfun(@(x) [x xs{ix}], v, 'UniformOutput', false);
                movefile(fns{:}, 'f');
            end
        end
        function pth   = ensureConsistentPwd(fqfp)
            %% ENSURECONSISTENTPWD
            %  @param fqfp is a filename or fileprefix
            %  @returns pth which is determined from fqfp; method cds to pth as needed
            
            pth = fileparts(fqfp);
            if (isempty(pth))
                pth = pwd;
                return
            end
            if (~strcmp(pwd, pth))
                cd(pth);
            end
        end
        function obj   = ensureLocalFourdfp(obj)
            %% ENSURELOCALFOURDFP
            %  @param obj is a data object to ensure to be in the pwd;
            %  supported:  cell, (Handle|)IOInterface, 4dfp, mgz, gz, nii.
            %  @return ensured contains the filename of desired data object in the pwd.
            
            this = mlfourdfp.FourdfpVisitor;            
            if (isempty(obj))
                obj = '';
                return
            end            
            if (iscell(obj))
                obj = cellfun(@(x) this.ensureLocalFourdfp(x), obj, 'UniformOutput', false);
                return
            end
            if (isa(obj, 'mlio.HandleIOInterface') || isa(obj, 'mlio.IOInterface'))
                obj.filepath = pwd;
                return
            end
            if (ischar(obj))
                obj_ = obj;
                [pth,fp,x] = myfileparts(obj);
                obj = [fp x];
                if (~this.isLocalFourdfp(obj_))
                    try
                        if (lexist_4dfp(fullfile(pth, fp)))
                            this.copyfile_4dfp(fullfile(pth, fp));
                            obj = [fp '.4dfp.ifh'];
                            return
                        end
                        if (strcmp(x, '.mgz') || strcmp(x, '.mgh'))
                            fp = ensureSafenameMgz(fp);
                            this.mri_convert(obj, [fp '.nii']);
                            this.nifti_4dfp_4(fp);
                            obj = [fp '.4dfp.ifh'];
                            return
                        end
                        if (lstrfind(x, '.gz'))
                            gunzip(obj, pwd);
                            this.nifti_4dfp_4(fp)
                            obj = [fp '.4dfp.ifh'];
                            return
                        end
                        if (strcmp(x, '.nii'))
                            this.nifti_4dfp_4(fullfile(pth, fp));
                            this.movefile_4dfp(fullfile(pth, fp));
                            obj = [fp '.4dfp.ifh'];
                            return
                        end
                    catch ME
                        dispexcept(ME);
                    end
                end
            end
        end
        function obj   = ensureSafeFileprefix(obj)
            %  @param obj may be empty; returns empty char.
            %  @param obj may be cell; if empty returns {''};
            %  @param obj may be mlio*IOInterface or char.
            
            this = mlfourdfp.FourdfpVisitor; 
            if (isempty(obj))
                obj = '';
                return
            end
            if (iscell(obj))
                obj = cellfun(@(x) this.ensureSafeFileprefix(x), obj, 'UniformOutput', false);
                return
            end
            if (isa(obj, 'mlio.IOInterface') || isa(obj, 'mlio.HandleIOInterface'))
                obj.fileprefix = this.ensureSafeFileprefix(obj.fileprefix);
                return
            end
            if (ischar(obj))
                [pth,fp,x] = myfileparts(obj);
                fn = [fp x];
                unsafe = this.UNSAFE_FILEPREFIXES;
                for u = 1:length(unsafe)
                    if (lstrfind(fn, unsafe{u}))
                        idxs = regexp(fn, unsafe{u});
                        while (~isempty(idxs))
                            fn = this.safesprintf(unsafe{u}, fn, idxs);
                            idxs = regexp(fn, unsafe{u});
                        end
                    end
                end
                obj0 = obj;
                obj = fullfile(pth, fn);
                if (lexist_4dfp(obj0) && ~lexist_4dfp(obj))
                    this.copyfile_4dfp(obj0, obj);
                end
                return
            end            
            error('mlfourdfp:unsupportedTypeclass', ...
                'FourdfpVisitor.ensureSafeFileprefix received %s', class(obj));
        end
        function m     = ifhMatrixSize(fqfn)
            if (~lstrfind(fqfn, '.4dfp.ifh'))
                fqfn = [fqfn '.4dfp.ifh'];
            end
            assert(lexist(fqfn));
            ifhp       = mlfourdfp.IfhParser.load(fqfn);
            m          = [0 0 0 0];
            [m(1),idx] = ifhp.rightSideNumeric('matrix size [1]');
            [m(2),idx] = ifhp.rightSideNumeric('matrix size [2]', idx);
            [m(3),idx] = ifhp.rightSideNumeric('matrix size [3]', idx);
             m(4)      = ifhp.rightSideNumeric('matrix size [4]', idx);
        end
        function tf    = isLocalFourdfp(obj)
            %  @param obj is a fourdfp data object:  filename, fileprefix, mlio.IOInterface, mlio.HandleIOInterface.
            %  @return tf is boolean for presence of fourdfp data object in the pwd.            

            if (isa(obj, 'mlio.HandleIOInterface') || isa(obj, 'mlio.IOInterface'))
                obj = obj.fqfilename;
            end            
            [~,fp] = myfileparts(obj);
            tf = lexist_4dfp(fp);
        end
        function tf    = lexist_4dfp(fqfp, varargin)
            %  @param fqfp exist as a filename on the filesystem or 
            %         fqfp is the fileprefix for 4dfp files on the filesystem.
            %  @return tf := true if useful files exist.
            
            tf = lexist(fqfp, varargin{:});
            if (~tf)
                tf = lexist([fqfp '.4dfp.hdr']) && ...
                     lexist([fqfp '.4dfp.ifh']) && ...
                     lexist([fqfp '.4dfp.img']);
            end
        end
        function tf    = lexist_mhdr(fqfp)
            tf = lexist(fqfp);
            if (~tf)
                [p,f] = myfileparts(fqfp);
                tf = lexist(fullfile(p, [f '.mhdr']));
            end
        end
        function tf    = lexist_vhdr(fqfp)
            tf = lexist(fqfp);
            if (~tf)
                [p,f] = myfileparts(fqfp);
                tf = lexist(fullfile(p, [f '.v.hdr']));
            end
        end
        function [s,r] = lns(fqfn)
            s = 0; r = '';
            if (~lexist(fullfile(pwd, mybasename(fqfn)), 'file') && ...
                 lexist(fqfn))
                [s,r] = mlbash(sprintf('ln -s %s', fqfn));
            end
            fprintf('mlfourdfp.FourdfpVisitor.lns->%s\n', fqfn);
        end
        function [s,r] = lns_4dfp(varargin)
            ip = inputParser;
            addRequired(ip, 'fqfpSrc',      @ischar);
            addOptional(ip, 'fqfpDest', '', @ischar);
            parse(ip, varargin{:});            
            fqfpSrc  = myfileprefix(ip.Results.fqfpSrc);
            fqfpDest = myfileprefix(ip.Results.fqfpDest);
            if (isempty(fqfpDest))
                fqfpDest = fullfile(pwd, mybasename(fqfpSrc));
            end     
            
            s = 0; r = '';
            ext = { '.hdr' '.ifh' '.img' '.img.rec' };         
            for e = 1:length(ext) 
                try
                    mlbash(sprintf('rm -f %s.4dfp%s', fqfpDest, ext{e}));
                catch ME
                    dispwarning(ME);
                end
                try
                    [s,r] = mlbash(sprintf('ln  -s %s.4dfp%s %s.4dfp%s', fqfpSrc, ext{e}, fqfpDest, ext{e}));
                catch ME
                    dispwarning(ME);
                end
            end
            fprintf('mlfourdfp.FourdfpVisitor.lns_4dfp->%s\n', fqfpDest);
        end
        function fn    = mri_convert(varargin)
            fn = mlsurfer.SurferVisitor.mri_convert(varargin{:});
        end
        function [s,r] = mkdir(pth)
            s = 0; r = '';
            if (~isdir(pth))
                [s,r] = mlbash(sprintf('mkdir -p %s', pth));
            end
        end
        function [s,r] = pushd(pth)
            s = 0; r = ''; %#ok<NASGU>
            [s,r] = mlbash(sprintf('pushd %s', pth));
        end
        function [s,r] = popd
            s = 0; r = ''; %#ok<NASGU>
            [s,r] = mlbash('popd');
        end
    end

	methods
        function [fqfp,s,r] = CT2mpr_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'mpr',              @this.lexist_4dfp);
            addRequired( ip, 'ct',               @this.lexist_4dfp);
            addParameter(ip, 'options', '',      @ischar);
            addParameter(ip, 'log', '/dev/null', @ischar);
            parse(ip, varargin{:});
            
            cmd = '%s %s %s';
            if (~strcmp(ip.Results.log, '/dev/null'))
                cmd = [cmd ' &> ' ip.Results.log]; 
            end
            [s,r] = this.CT2mpr_4dfp__(sprintf( ...
                cmd, ip.Results.mpr, ip.Results.ct, ip.Results.options));
            fqfp = fullfile(pwd, sprintf('%s_on_%s', mybasename(ip.Results.ct), mybasename(ip.Results.mpr)));
        end        
        function      [s,r] = C2S_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in', @this.lexist_4dfp);
            addOptional( ip, 'out', [varargin{1} 'S'], @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.C2T_4dfp__( ...
                sprintf('%s %s', ip.Results.in, [ip.Results.in 'T'])); %#ok<ASGLU>
            [s,r] = this.T2S_4dfp__( ...
                sprintf('%s %s', [ip.Results.in 'T'], ip.Results.out));
        end
        function      [s,r] = C2T_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',      @this.lexist_4dfp);
            addOptional( ip, 'out', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.C2T_4dfp__( ...
                sprintf('%s %s', ip.Results.in, ip.Results.out));
        end
        function      [s,r] = IFhdr_to_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'IFstr',                              @this.lexist_vhdr);
            addOptional(ip, 'outroot', myfileprefix(varargin{1}), @ischar);
            parse(ip, varargin{:});
                    
            [p,f] = myfileparts(ip.Results.IFstr);
            [s,r] = this.IFhdr_to_4dfp__(sprintf( ...
                '%s %s', fullfile(p, [f '.v.hdr']), ip.Results.outroot));   
            
            delete([ip.Results.outroot 'fz*']);
        end
        function      [s,r] = S2C_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in', @this.lexist_4dfp);
            addOptional( ip, 'out', [varargin{1} 'C'], @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.S2T_4dfp__( ...
                sprintf('%s %s', ip.Results.in, [ip.Results.in 'T'])); %#ok<ASGLU>
            [s,r] = this.T2C_4dfp__( ...
                sprintf('%s %s', [ip.Results.in 'T'], ip.Results.out));
        end
        function      [s,r] = S2T_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',      @this.lexist_4dfp);
            addOptional( ip, 'out', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.S2T_4dfp__( ...
                sprintf('%s %s', ip.Results.in, ip.Results.out));
        end
        function      [s,r] = T2C_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',      @this.lexist_4dfp);
            addOptional( ip, 'out', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.T2C_4dfp__( ...
                sprintf('%s %s', ip.Results.in, ip.Results.out));
        end
        function      [s,r] = T2S_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',      @this.lexist_4dfp);
            addOptional( ip, 'out', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.T2S_4dfp__( ...
                sprintf('%s %s', ip.Results.in, ip.Results.out));
        end
        function      [s,r] = actmapf_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'format',      @ischar);
            addRequired( ip, 'input',       @this.lexist_4dfp);
            addParameter(ip, 'options', '', @ischar);
            addParameter(ip, 'normalized', false, @islogical);
            parse(ip, varargin{:});
            
            if (ip.Results.normalized)
                this.normalizeFrames(ip.Results.input);
            end            
            [s,r] = this.actmapf_4dfp__(sprintf( ...
                '%s %s %s', ip.Results.format, ip.Results.input, ip.Results.options));
        end        
        function [t4,fqfp,s,r] = ...
                              align_2051(this, varargin)
            %% ALIGN_2051 calls imgreg_4dfp with modes [2051 2051] and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_modes2, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            [t4,fqfp,s,r]  = this.align_modes2(2051, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_3075(this, varargin)
            %% ALIGN_3075 calls imgreg_4dfp with mode 3075 and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_mode, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            [t4,fqfp,s,r]  = this.align_mode(3075, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_4099(this, varargin)
            %% ALIGN_4099 calls imgreg_4dfp with modes [4099 4099] and writes a log.
            
            [t4,fqfp,s,r]  = this.align_modes2(4099, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_10243(this, varargin)
            %% ALIGN_10243 calls imgreg_4dfp with modes [10243 10243] and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_modes2, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            [t4,fqfp,s,r]  = this.align_modes2(10243, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_commonModal(this, varargin)
            %% ALIGN_COMMONMODAL calls imgreg_4dfp with modes [4099 4099 3075 2051 2051 10243] + 256 and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_modes6, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            %assert(~lstrfind(varargin, 'useCommonModal'), ...
            %    'mlfourdfp:ambiguousVargin', ...
            %    'FourdfpVisitor.align_multiSpectral.varagin.useCommonModal:multiplyAssigned');
            [t4,fqfp,s,r] = this.align_modes6('useCommonModal', true, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_commonModal7(this, varargin)
            %% ALIGN_COMMONMODAL7 calls imgreg_4dfp with modes [4099 4099 3075 2051 2051 10243] + 256 and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_modes7, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            %assert(~lstrfind(varargin, 'useCommonModal'), ...
            %    'mlfourdfp:ambiguousVargin', ...
            %    'FourdfpVisitor.align_multiSpectral.varagin.useCommonModal:multiplyAssigned');
            [t4,fqfp,s,r] = this.align_modes7('useCommonModal', true, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_crossModal(this, varargin)
            %% ALIGN_CROSSMODAL calls imgreg_4dfp with modes [4099 4099 3075 2051 2051 10243] and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_modes6, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            %assert(~lstrfind(varargin, 'useCommonModal'), ...
            %    'mlfourdfp:ambiguousVargin', ...
            %    'FourdfpVisitor.align_multiSpectral.varagin.useCommonModal:multiplyAssigned');
            [t4,fqfp,s,r] = this.align_modes6('useCommonModal', false, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_crossModal7(this, varargin)
            %% ALIGN_CROSSMODAL7 calls imgreg_4dfp with modes [4099 4099 3075 2051 2051 10243] and writes a log.
            %  See also mlfourdfp.FourdfpVisitor.align_modes7, mlfourdfp.FourdfpVisitor.imgreg_4dfp.
            
            %assert(~lstrfind(varargin, 'useCommonModal'), ...
            %    'mlfourdfp:ambiguousVargin', ...
            %    'FourdfpVisitor.align_multiSpectral.varagin.useCommonModal:multiplyAssigned');
            [t4,fqfp,s,r] = this.align_modes7('useCommonModal', false, varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_multiSpectral(this, varargin)
            [t4,fqfp,s,r] = this.align_crossModal(varargin{:});
        end
        function [t4,fqfp,s,r] = ...
                              align_TOF(this, varargin)
            %% ALIGN_TOF calls imgreg_4dfp with modes [4099 4099 3075 2051 2051 10243] and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param destMask   "
            %  @param source     "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param useMetricGradient:  is logical; cf. ${TRANSFER}/cross_modal_intro.pdf.
            %  @param t4img_4dfp is logical.
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            ip = inputParser;
            addParameter(ip, 'dest',       '',     @this.lexist_4dfp);
            addParameter(ip, 'source',     '',     @this.lexist_4dfp);
            addParameter(ip, 'destMask',   'none', @ischar);
            addParameter(ip, 'sourceMask', 'none', @ischar);
            addParameter(ip, 'destBlur',   0,      @isnumeric);
            addParameter(ip, 'sourceBlur', 0,      @isnumeric);
            addParameter(ip, 't40',        this.transverse_t4, @(x) lexist(x, 'file')); 
            addParameter(ip, 't4',         '',     @ischar);
            addParameter(ip, 'useMetricGradient', false, @islogical);
            addParameter(ip, 't4img_4dfp',        true, @islogical);
            addParameter(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});
            dest       = ip.Results.dest;
            source     = ip.Results.source;
            destMask   = ip.Results.destMask;
            sourceMask = ip.Results.sourceMask;
            t4         = ip.Results.t4;
            log        = ip.Results.log;
            
            if (strcmp(dest, source))
                t4 = this.filenameT4(source, dest);
                copyfile(this.transverse_t4, t4); % identity transformation
                fqfp = '';
                return
            end
            if (isempty(t4))
                t4 = this.filenameT4(source, dest); end
            if (ip.Results.destBlur > 0)
                dest = this.imgblur_4dfp(dest, ip.Results.destBlur); end
            if (ip.Results.sourceBlur > 0)
                source = this.imgblur_4dfp(source, ip.Results.sourceBlur); end
            
            if (~strcmp(ip.Results.t40, t4))
                [s,r] = copyfile(ip.Results.t40, t4, 'f'); %#ok<ASGLU>
            end
            [s,r] = dbbash(sprintf('chmod 777 %s', t4)); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 4099, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 4099, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 3075, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051, log); %#ok<ASGLU>       
            if (ip.Results.t4img_4dfp)
                [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'options', ['-O' ip.Results.dest]);
            else
                fqfp = ''; s = 0; r = '';
            end
        end
        function      [s,r] = compute_defined_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'in',          @this.lexist_4dfp);
            addOptional(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.compute_defined_4dfp__( ...
                sprintf(' %s %s', ip.Results.options, ip.Results.in));
        end
        function ic         = convertImageToLocal4dfp(this, varargin)
            %% CONVERTIMAGETOLOCAL4DFP
            %  @param required fqfn is char.
            %  @param optional fqfn1 is char and may be lacking a file extension.
            %  @return ic is an mlfourd.ImagingContext for 4dfp data created in the the cwd.
            
            ip = inputParser;
            addRequired(ip, 'fqfn', @(x) lexist(x, 'file'));
            addOptional(ip, 'fqfn1', '', @ischar);
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            if (isempty(ipr.fqfn1))
                ic = this.convertImageToLocal4dfp_1arg(ipr);
            else
                ic = this.convertImageToLocal4dfp_2arg2(ipr);
            end
        end
        function ic         = convertImageToLocal4dfp_1arg(this, ipr)
            % See also:  convertImageToLocal4dfp
            assert(isstruct(ipr));
            [~,fp] = myfileparts(ipr.fqfn);
            locfp  = fullfile(pwd, fp);
            nii    = [locfp '.nii'];
            ifh    = [locfp '.4dfp.ifh'];
            if (~lexist_4dfp(ifh))
                mlsurfer.SurferVisitor.mri_convert(ipr.fqfn, nii);
                this.nifti_4dfp_4(locfp);
            end
            ic = mlfourd.ImagingContext(ifh);
        end
        function ic         = convertImageToLocal4dfp_2arg2(this, ipr)
            % See also:  convertImageToLocal4dfp
            assert(isstruct(ipr));
            [~,fp]  = myfileparts(ipr.fqfn);
            [~,fp1] = myfileparts(ipr.fqfn1);
            locfp   = fullfile(pwd, fp);
            locfp1  = fullfile(pwd, fp1);
            nii     = [locfp '.nii'];
            ifh     = [locfp '.4dfp.ifh'];
            ifh1    = [locfp1 '.4dfp.ifh'];
            if (~lexist_4dfp(ifh1))
                mlsurfer.SurferVisitor.mri_convert(ipr.fqfn, nii);
                this.nifti_4dfp_4(locfp);
                this.move_4dfp(locfp, locfp1);
            end
            ic = mlfourd.ImagingContext(ifh);
        end        
        function      [s,r] = copy_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'in',  @ischar);
            addRequired(ip, 'out', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.copy_4dfp__( ...
                sprintf(' %s %s', ip.Results.in, ip.Results.out));
        end
        function      [s,r] = cropfrac_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'frac', @isnumeric);
            addRequired( ip, 'in',   @this.lexist_4dfp);
            addRequired( ip, 'out',  @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.cropfrac_4dfp__( ...
                sprintf(' %g %s %s', ip.Results.frac, ip.Results.in, ip.Results.out));
        end
        function      [s,r] = dcm_sort(this, varargin)
            ip = inputParser;
            addRequired( ip, 'pth',         @isdir);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.dcm_sort__( ...
                sprintf(' %s %s', ip.Results.pth, ip.Results.options));
        end
        function      [s,r] = dcm_sort_PPG(this, varargin)
            ip = inputParser;
            addRequired( ip, 'pth',         @isdir);
            parse(ip, varargin{:});
            
            [s,r] = this.dcm_sort_PPG__( ...
                sprintf(' %s', ip.Results.pth));
        end
        function      [s,r] = dcm_to_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'target',          @(x) isdir(fileparts(x))); %#ok<*ISDIR>
            addParameter(ip, 'base', 'analyze', @ischar);
            addParameter(ip, 'options', '',     @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.dcm_to_4dfp__( ...
                sprintf(' -b %s %s %s', ip.Results.base, ip.Results.options, ip.Results.target));
        end
        function      [s,r] = delete_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',  @this.lexist_4dfp);
            parse(ip, varargin{:});
            
            [s,r] = this.delete_4dfp__( ...
                sprintf(' %s', ip.Results.in));
        end
        function      [e,c] = etaAndCurves(this, varargin)
            ip = inputParser;
            addRequired(ip, 'dest',   @this.lexist_4dfp);
            addRequired(ip, 'source', @this.lexist_4dfp);
            addOptional(ip, 't4',     ['_' this.filenameT4(varargin{1}, varargin{2})], @ischar);
            addOptional(ip, 'mode',   2051, @isnumeric);
            addOptional(ip, 'log',    loggerFilename('', 'func', 'FourdfpVisitor_eta'), @ischar);
            parse(ip, varargin{:});
            
            this.imgreg_4dfp__(sprintf('%s none %s none %s %i >> %s', ...
                this.ensureImgregFiles(ip.Results.dest), ...
                this.ensureImgregFiles(ip.Results.source), ...
                ip.Results.t4, ...
                ip.Results.mode, ...
                ip.Results.log));
            lp = mlio.LogParser.load(ip.Results.log);
            c  = lp.nextLineNNumeric('100000*second partial in parameter space', 6, 1);
            e  = lp.rightSideNumeric('eta,q', lp.length);
            delete(ip.Results.t4);
            delete(ip.Results.log);
        end
        function      [s,r] = extract_frame_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp',        @(x) this.lexist_4dfp(x));
            addRequired(ip, 'frame',       @isnumeric);
            addOptional(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.extract_frame_4dfp__( ...
                sprintf('%s %i %s', ip.Results.fqfp, ip.Results.frame, ip.Results.options));
        end
        function    fn      = filenameT4(this, fn1, fn2)
            %% FILENAMET4
            %  @param fn1 may be fully-qualified; a filename or fileprefix; a 4dfp file or t4-file.
            %  @param fn2 "
            %  @returns fn is a t4-filename for transformations of fn1 to fn2 in the folder of fn2.
            %  N.B.:  fn1 and fn2 must both be 4dfp or both be t4.
            
            if (length(fn1) > 2 && length(fn2) > 2)
                if (strcmp(fn1(end-2:end), '_t4') && strcmp(fn2(end-2:end), '_t4'))
                    fn = this.filenameT4mul(fn1, fn2);
                    return
                end            
                if (strcmp(fn1(end-2:end), '_t4') || strcmp(fn2(end-2:end), '_t4'))
                    error('mlfourdfp:filenameError', 'FourdfpVisitor.filenameT4:  expected 4dfp but received t4');
                end
                if (lstrfind(fn1, '_to_') || lstrfind(fn2, '_to_'))
                    error('mlfourdfp:filenameError', 'FourdfpVisitor.filenameT4:  malformed filename:  %s, %s', fn1, fn2);
                end
            end           
            
            [~, f1] = myfileparts(fn1);
            [p2,f2] = myfileparts(fn2);
            [first,last] = this.parseFirstOnLastBasenames(f1, f2);
            fn = fullfile(p2, sprintf('%s_to_%s_t4', first, last));
        end
        function    fn      = filenameT4inv(~, fn)
            %% FILENAMET4INV
            %  @param fn is a t4-filename 
            %  @returns a t4-filename appropriate for the inverse transformation
            
            [p,f] = myfileparts(fn);
            if (~lstrfind(fn, '_t4'))
                error('mlfourdfp:filenameError', 'FourdfpVisitor.filenameT4inv received a malformed t4 filename:  %s', fn);
            end
            r  = regexp(f, '(?<first>\S+)_to_(?<last>\S+)_t4', 'names');
            fn = fullfile(p, sprintf('%s_to_%s_t4', r.last, r.first));
        end
        function    fn      = filenameT4mul(this, fn1, fn2)
            %% FILENAMET4MUL
            %  @param fn1 may be fully-qualified; is a t4-filename.
            %  @param fn2 "
            %  @returns fn is a t4-filename for composite transformations of fn1 to fn2 in the folder of fn2.
            
            [~, f1] = myfileparts(fn1);
            [p2,f2] = myfileparts(fn2);
            if (strcmp(f1, 'T_t4'))
                fn = fullfile(p2, f2);
                return
            end
            if (strcmp(f2, 'T_t4'))
                fn = fullfile(p2, f1);
                return
            end
            [first,last] = this.parseFirstToLastBasenames(f1, f2);
            fn = fullfile(p2, sprintf('%s_to_%s_t4', first, last));
        end
        function  fqfp      = fileprefixBlurred(~, fqfp, blurArg)
            if (prod(blurArg) < eps)
                return
            end
            fqfp = sprintf('%s_b%i', fqfp, round(blurArg*10));
        end
        function  fqfp      = fileprefixGaussed(~, fqfp, gaussArg)
            fqfp = sprintf('%s_g%i', fqfp, round(gaussArg*10));
        end
        function  fqfp      = fileprefixT4img(this, fp1, fp2)
            %% FILENAMET4IMG
            %  @param fp1 may be fully-qualified; a 4dfp fileprefix or filename or t4-file.
            %  @param fp2 "
            %  @returns fqfp is a 4dfp-filename for fp1 on fp2 in the folder of fn2.
            
            [p1,f1] = myfileparts(fp1);
            [~, f2] = myfileparts(fp2);
            [first,last] = this.parseFirstOnLastBasenames(f2, f1);
            fqfp = fullfile(p1, sprintf('%s_on_%s', first, last));
        end
        function [fqfp,s,r] = flip_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'basis',   @ischar);
            addRequired(ip, 'input',   @this.lexist_4dfp);
            addOptional(ip, 'outroot', sprintf('%s_flip%s', varargin{2}, varargin{1}), @ischar);
            parse(ip, varargin{:});
            
            %if (~(strcmp('x', ip.Results.basis) || strcmp('y', ip.Results.basis) || strcmp('z', ip.Results.basis)))
            %    return
            %end
            fqfp = ip.Results.outroot;            
            [s,r] = this.flip_4dfp__(sprintf('-%s %s %s', ip.Results.basis, ip.Results.input, ip.Results.outroot));
        end
        function      [s,r] = freesurfer2mpr_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'mpr',              @ischar);
            addRequired( ip, 'fs',               @ischar);
            addParameter(ip, 'options', '',      @ischar);
            addParameter(ip, 'log', '/dev/null', @ischar);
            parse(ip, varargin{:});
            mpr_ = ip.Results.mpr;
            fs_  = ip.Results.fs;
            
            cmd = '%s %s %s setecho';
            if (~strcmp(ip.Results.log, '/dev/null'))
                cmd = [cmd ' &> ' ip.Results.log]; 
            end
            [s,r] = this.freesurfer2mpr_4dfp__(sprintf( ...
                cmd, mpr_, fs_, ip.Results.options));
        end
        function [fqfp,s,r] = gauss_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'input',   @this.lexist_4dfp);
            addRequired(ip, 'f_half',  @isnumeric);
            addOptional(ip, 'outroot', sprintf('%s_g%s', varargin{1}, num2str(varargin{2})), @ischar);
            parse(ip, varargin{:});
            
            if (ip.Results.f_half < eps)
                fqfp = ip.Results.input;
                s = 0; r = '';
                return
            end
            fqfp = ip.Results.outroot;            
            [s,r] = this.gauss_4dfp__(sprintf('%s %g %s', ip.Results.input, ip.Results.f_half, ip.Results.outroot));
        end
        function [fqfp,s,r] = imgblur_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp',        @this.lexist_4dfp);
            addRequired(ip, 'fwhm',        @isnumeric);
            addOptional(ip, 'options', '', @ischar);            
            parse(ip, varargin{:});

            if (ip.Results.fwhm < eps)
                fqfp = ip.Results.fqfp;
                s = 0; r = '';
                return
            end
            fqfp  = this.fileprefixBlurred(ip.Results.fqfp, ip.Results.fwhm);
            [s,r] = this.imgblur_4dfp__( ...
                sprintf('%s %s %g', ip.Results.options, ip.Results.fqfp, ip.Results.fwhm));
        end
        function      [s,r] = imgblur_all_4dfp(this, varargin)
            %% IMGBLUR_ALL_4DFP applies imgblur_4dfp to all 4dfp sets in the current working directory      

            eSet = mlsystem.DirTool('*.4dfp.ifh');
            for iSet = 1:length(eSet.fns)                                
                [s,r] = this.imgblur_4dfp(eSet.fns{iSet}, varargin{:});
            end
        end
        function      [s,r] = imgopr_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'operation',   @ischar);
            addRequired( ip, 'out',         @ischar);
            addRequired( ip, 'in1',         @this.lexist_4dfp);
            addRequired( ip, 'in2',         @this.lexist_4dfp);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.imgopr_4dfp__( ...
                sprintf('-%s%s %s %s %s', ip.Results.operation, ip.Results.out, ip.Results.in1, ip.Results.in2, ip.Results.options));
        end
        function      [s,r] = imgreg_4dfp(this, varargin)
            %% IMGREG_4DFP 
            %         mode word bit assignments
            %     1   enable coordinate transform
            %     2   enable 3D alignment
            %     4   enable affine warp (12 parameters in 3D 6 parameters in 2D)
            %     8   enable voxel size adjust
            %    16   disable x voxel size adjust
            %    32   disable y voxel size adjust
            %    64   disable z voxel size adjust
            %   128   unassigned
            %   256   when set use difference image minimization (common mode for similar contrast mechanisms)
            %   512   superfine mode (2 mm cubic grid metric sampling)
            %  1024   fast mode (12 mm cubic grid metric sampling)
            %  2048   fine mode (5 mm cubic grid metric sampling)
            %  4096   [T] restricted to translation explored at 7.5 mm intervals
            %  8192   enable parameter optimization by computation of the metric gradient
            %         in parameter space and inversion of the Hessian
            %
            %  @param dest is 4dfp, fileprefix or filename 
            %  @param destMask is 4dfp or 'none'
            %  @param source is 4dfp
            %  @param sourceMask is 4dfp or 'none'
            %  @param t4 exists on the filesystem
            %  @param mode is numeric (default := 2051)
            %  @param log is on the filesystem (default := '/dev/null')

            ip = inputParser;
            addRequired(ip, 'dest',       @this.lexist_4dfp);
            addRequired(ip, 'destMask',   @ischar); % normally unnecessary
            addRequired(ip, 'source',     @this.lexist_4dfp);
            addRequired(ip, 'sourceMask', @ischar);
            addOptional(ip, 't4',         this.filenameT4(varargin{1}, varargin{3}), @ischar);
            addOptional(ip, 'mode',       2051, @isnumeric);
            addOptional(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.imgreg_4dfp__(sprintf('%s %s %s %s %s %i >> %s', ...
                this.ensureImgregFiles(ip.Results.dest), ...
                this.ensureImgregFiles(ip.Results.destMask), ...
                this.ensureImgregFiles(ip.Results.source), ...
                this.ensureImgregFiles(ip.Results.sourceMask), ...
                ip.Results.t4, ...
                ip.Results.mode, ...
                ip.Results.log));
        end
        function [fqfp,s,r] = maskimg_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'imgfile', @this.lexist_4dfp);
            addRequired(ip, 'mskfile', @this.lexist_4dfp);
            addRequired(ip, 'outfile', @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.maskimg_4dfp__(sprintf('%s %s %s %s', ...
                ip.Results.options, ip.Results.imgfile, ip.Results.mskfile, ip.Results.outfile));
            fqfp = ip.Results.outfile;
        end
        function      [s,r] = move_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',  @this.lexist_4dfp);
            addRequired( ip, 'out', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.copy_4dfp__( ...
                sprintf(' %s %s', ip.Results.in, ip.Results.out));
            x = {'ifh' 'hdr' 'img' 'img.rec'};
            for ix = 1:length(x)
                assert(lexist( sprintf(   '%s.4dfp.%s', ip.Results.out, x{ix})));
                [s,r] = mlbash(sprintf('rm %s.4dfp.%s', ip.Results.in,  x{ix}));
            end
        end
        function      [s,r] = mpr2atl_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',               @ischar);
            addParameter(ip, 'options', sprintf('-T%s/TRIO_Y_NDC -S711-2B', getenv('REFDIR')), @ischar);
            addParameter(ip, 'log', '/dev/null', @ischar);
            parse(ip, varargin{:});
            in_ = ip.Results.in;
            
            cmd = '%s %s';
            if (~strcmp(ip.Results.log, '/dev/null'))
                cmd = [cmd ' &> ' ip.Results.log]; 
            end
            [s,r] = this.mpr2atl_4dfp__(sprintf( ...
                cmd, in_, ip.Results.options));
        end
        function      [s,r] = msktgen_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',             @ischar);
            addOptional( ip, 'threshold', 200, @isnumeric);
            addParameter(ip, 'options',  '',   @ischar);
            addParameter(ip, 'log',      '',   @ischar);
            parse(ip, varargin{:});
            log = ip.Results.log;
            if (~isempty(ip.Results.log))
                log = [' &> ' log];
            end           
            in_ = ip.Results.in;
            [s,r] = this.msktgen_4dfp__( ...
                sprintf('%s %i %s %s', in_, ip.Results.threshold, ip.Results.options, log));
        end
        function      [s,r] = msktgen_b110_4dfp(this, varargin)
            %% MSKTGEN_B110_4DFP uses 11 mm fwhh blur of a mask of the whole brain.
            
            ip = inputParser;
            addRequired( ip, 'in',           @ischar);
            addOptional( ip, 'threshold', 0, @isnumeric);
            addParameter(ip, 'options',  '', @ischar);
            addParameter(ip, 'log',      '', @ischar);
            parse(ip, varargin{:});
            log = ip.Results.log;
            if (~isempty(ip.Results.log))
                log = [' &> ' log];
            end            
            in_ = ip.Results.in;
            [s,r] = this.msktgen_b110_4dfp__( ...
                sprintf('%s %i %s %s', in_, ip.Results.threshold, ip.Results.options, log));
        end
        function      [s,r] = msktgen2_4dfp(this, varargin)
            %% MSKTGEN2_4DFP masks the entire head for TRIO_Y_NDC.
            
            ip = inputParser;
            addRequired( ip, 'in',           @ischar);
            addOptional( ip, 'threshold', 0, @isnumeric);
            addParameter(ip, 'options',  '', @ischar);
            addParameter(ip, 'log',      '', @ischar);
            parse(ip, varargin{:});
            log = ip.Results.log;
            if (~isempty(ip.Results.log))
                log = [' &> ' log];
            end            
            in_ = ip.Results.in;
            [s,r] = this.msktgen2_4dfp__( ...
                sprintf('%s %i %s %s', in_, ip.Results.threshold, ip.Results.options, log));
        end
        function      [s,r] = msktgen3_4dfp(this, varargin)
            %% MSKTGEN3_4DFP masks in the pharynx for TRIO_Y_NDC.
            
            ip = inputParser;
            addRequired( ip, 'in',             @ischar);
            addOptional( ip, 'threshold', 200, @isnumeric);
            addParameter(ip, 'options',   '',  @ischar);
            addParameter(ip, 'log',       '',  @ischar);
            parse(ip, varargin{:});
            log = ip.Results.log;
            if (~isempty(ip.Results.log))
                log = [' &> ' log];
            end
            in_ = ip.Results.in;          
            [s,r] = this.msktgen3_4dfp__( ...
                sprintf('%s %i %s %s', in_, ip.Results.threshold, ip.Results.options, log));
        end
        function this       = msktgenMprage(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfp', @ischar);
            addOptional(ip, 'atl', fullfile(getenv('REFDIR'), 'TRIO_Y_NDC'), @ischar);
            parse(ip, varargin{:});
            fqfp = ip.Results.fqfp;
            atl  = ip.Results.atl;
            
            log = sprintf('msktgenMprage_%s.log', datestr(now, 30));
            this.mpr2atl_4dfp(fqfp, 'options', ['-T' atl], 'log', log);
            this.msktgen_4dfp(fqfp, 'options', ['-T' atl], 'log', log);
        end   
        function      [s,r] = nifti_4dfp_4(this, varargin)
            %% NIFTI_4DFP_4 converts .nii[.gz] to .4dfp.*, deleting the .nii[.gz] afterwards.
            %  @param fileprefix is char.
            %  @param fileprefixOut is char.
            %  @param named minusN is logical; true =: sends -N flag to nifti_4dfp to remove center parameters from
            %  *.4dfp.ifh.
            
            ip = inputParser;
            addRequired(ip, 'fileprefix', @ischar);
            addOptional(ip, 'fileprefixOut', varargin{1}, @ischar);
            addParameter(ip, 'minusN', true);
            parse(ip, varargin{:}); 
            [pth,fp] = myfileparts(ip.Results.fileprefix);
            fp = fullfile(pth, fp);
            [ptho,fpo] = myfileparts(ip.Results.fileprefixOut);
            fpo = fullfile(ptho, fpo);
            
            % manage .nii.gz
            if (lexist([fp '.nii.gz'], 'file') && ...
               ~lexist([fp '.nii'],    'file'))
                gunzip([fp '.nii.gz']);
            end
            if (lexist([fp '.nii.gz'], 'file') && ...
                lexist([fp '.nii'],    'file'))
                delete([fp '.nii.gz']);
            end
            
            if (lexist([fp '.nii']))
                if (~ip.Results.minusN && ...
                    (lstrfind(fp, '111') || lstrfind(fp, '222') || lstrfind(fp, '333') || ...
                     lstrfind(fp, 'TRIO_Y_NDC') || ...
                     lstrfind(fp, '711-2') || ...
                     lstrfind(fp, '_atlas')))
                    [s,r] = this.nifti_4dfp__(sprintf(' -4 %s.nii %s.4dfp.ifh', fp, fpo));
                    deleteExisting([fp '.nii']);
                    return
                end
                [s,r] = this.nifti_4dfp__(sprintf(' -4 %s.nii %s.4dfp.ifh -N', fp, fpo));
                deleteExisting([fp '.nii']);
                return
                
                %deleteExisting([fp '.4dfp.img_to_atlas_t4']); % incipient BUG
            end
            error('mlfourdfp:fileNotFound', 'FourdfpVisitor.nifti_4dfp_4:  %s.nii not found', fp);
        end
        function      [s,r] = nifti_4dfp_n(this, varargin)
            %% NIFTI_4DFP_4 converts .4dfp.* to .nii[.gz], keeping all .4dfp.* afterwards.
            %  @param fileprefix is char.
            %  @param fileprefixOut is char.
            
            ip = inputParser;
            addRequired(ip, 'fileprefix',                 @ischar);
            addOptional(ip, 'fileprefixOut', varargin{1}, @ischar);
            parse(ip, varargin{:});
            [pth,fp] = myfileparts(ip.Results.fileprefix);
            fp = fullfile(pth, fp);
            [ptho,fpo] = myfileparts(ip.Results.fileprefixOut);
            fpo = fullfile(ptho, fpo);            

            if (lexist([fp '.4dfp.ifh']))
                [s,r] = this.nifti_4dfp__(sprintf(' -n %s.4dfp.ifh %s.nii', fp, fpo));
                return
            end
            error('mlfourdfp:fileNotFound', 'FourdfpVisitor.nifti_4dfp_n:  %s.4dfp.* not found', fp);
        end
        function      [s,r] = nifti_4dfp_ng(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fileprefix',                 @ischar);
            addOptional(ip, 'fileprefixOut', varargin{1}, @ischar);
            parse(ip, varargin{:});
            [pth,fp] = myfileparts(ip.Results.fileprefix);
            fp = fullfile(pth, fp);
            [ptho,fpo] = myfileparts(ip.Results.fileprefixOut);
            fpo = fullfile(ptho, fpo);            

            if (lexist([fp '.4dfp.ifh']))
                [s,r] = this.nifti_4dfp__(sprintf(' -n %s.4dfp.ifh %s.nii', fp, fpo));
                gzipExisting(  [fpo '.nii']);
                deleteExisting([fpo '.nii']);
                return
            end
            error('mlfourdfp:fileNotFound', 'FourdfpVisitor.nifti_4dfp_ng:  %s.4dfp.* not found', fp);
        end        
        function              normalizeFrames(this, fqfp)
            tmp = tmpFileprefix('func', 'normalizeFrames');
            this.move_4dfp(fqfp, tmp);            
            
            this.nifti_4dfp_ng(tmp);
            ic = mlfourd.ImagingContext(tmp);
            nii = ic.niftid;
            nii.img = this.normalizeFramesBySums(nii.img);            
            ic = mlfourd.ImagingContext(nii);
            ic.saveas(fqfp);
            this.nifti_4dfp_4(fqfp);
            delete([fqfp '.nii*']);
            delete([fqfp '.log']);
            delete([tmp  '.nii*']);
            delete([tmp  '.4dfp.*']);
        end
        function      [f,l] = parseFilenameT4(~, fn)
            [pth,fp,x] = fileparts(fn);
            assert(isempty(x), 'unexpected filesuffix %s', x);
            
            if (~contains(fp, '_t4'))
                error('mlfourdfp:filenameError', 'FourdfpVisitor.parseFilenameT4 received a malformed t4 filename:  %s', fn);
            end
            r = regexp(fp, '(?<first>\S+)_to_(?<last>\S+)_t4', 'names');
            f = fullfile(pth, r.first);
            l = r.last;
        end
        function      [f,l] = parseFirstLastBasenames(~, name1, name2, sep)
            name1 = mybasename(name1);
            name2 = mybasename(name2);    
            assert(ischar(sep));
            if (lstrfind(name1, sep))
                idx = strfind(name1, sep);
                f   = name1(1:idx(1)-1);
            elseif (contains(name1, '_to_'))
                idx = strfind(name1, '_to_');
                f   = name1(1:idx(1)-1);
            else
                f   = name1;
            end
            if (lstrfind(name2, sep))
                idx = strfind(name2, sep);
                l   = name2(idx(end)+4:end);
            elseif (lstrfind(name2, '_to_'))
                idx  = strfind(name2, '_to_');
                idx2 = strfind(name2, '_t4');
                l    = name2(idx(end)+4:idx2(end)-1);
            else
                l   = name2;
            end            
        end
        function      [f,l] = parseFirstOnLastBasenames(this, name1, name2)
            [f,l] = this.parseFirstLastBasenames(name1, name2, '_on_');
        end
        function      [f,l] = parseFirstOpLastBasenames(this, name1, name2)
            [f,l] = this.parseFirstLastBasenames(name1, name2, '_op_');
        end
        function      [f,l] = parseFirstToLastBasenames(~, name1, name2)
            name1 = mybasename(name1);
            name2 = mybasename(name2);
            if (lstrfind(name1, '_to_'))
                idx = strfind(name1, '_to_');
                f   = name1(1:idx(1)-1);
            else
                error('mlfourdfp:filenameError',  ...
                      'FourdfpVisitor.parseFirstToLastBasenames received malformed name->%s', name1);
            end
            if (lstrfind(name2, '_to_'))
                idx  = strfind(name2, '_to_');
                idx2 = strfind(name2, '_t4');
                l    = name2(idx(end)+4:idx2(end)-1);
            else
                error('mlfourdfp:filenameError',  ...
                      'FourdfpVisitor.parseFirstToLastBasenames received malformed name->%s', name2);
            end            
        end
        function      [s,r] = paste_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'inlist',      @lexist);
            addRequired( ip, 'outfile',     @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.paste_4dfp__( ...
                sprintf(' %s %s %s', ip.Results.options, ip.Results.inlist, ip.Results.outfile));
        end
        function      [s,r] = pseudo_dcm_sort(this, varargin)
            ip = inputParser;
            addRequired( ip, 'path',        @isdir);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.pseudo_dcm_sort__( ...
                sprintf(' %s %s', ip.Results.path, ip.Results.options));
        end
        function      [s,r] = scale_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',          @this.lexist_4dfp);
            addRequired( ip, 'scale',       @isnumeric);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.scale_4dfp__( ...
                sprintf('%s %g %s', ip.Results.in, ip.Results.scale, ip.Results.options));
        end
        function [fqfp,s,r] = sqrt_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',          @this.lexist_4dfp);
            addOptional( ip, 'outroot', '', @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.sqrt_4dfp__( ...
                sprintf('%s %s %s', ip.Results.in, ip.Results.outroot, ip.Results.options));
            if (~isempty(ip.Results.outroot))
                fqfp = ip.Results.outroot;
            else
                fqfp = [ip.Results.in '_sqrt'];
            end
        end
        function      [s,r] = sif_4dfp(this, varargin)
            ip = inputParser;
            addRequired(ip, 'sifstr',      @this.lexist_mhdr);
            addOptional(ip, 'outroot', '', @ischar);
            parse(ip, varargin{:});
            outroot = ip.Results.outroot;
            
            [p,f] = myfileparts(ip.Results.sifstr);
            if (isempty(outroot)); outroot = fullfile(p, f); end
            [s,r] = this.sif_4dfp__(sprintf( ...
                '%s %s', fullfile(p, [f '.mhdr']), outroot));
        end
        function [fqfp,s,r] = t4img_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 't4',          @lexist);
            addRequired( ip, 'in',          @this.lexist_4dfp);
            addParameter(ip, 'out', this.fileprefixT4img(varargin{1:2}), ...
                                            @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.t4img_4dfp__( ...
                sprintf('%s %s %s %s', ...
                ip.Results.t4, ...
                myfileprefix(ip.Results.in), ...
                myfileprefix(ip.Results.out), ...
                ip.Results.options));
            fqfp = myfileprefix(ip.Results.out);
        end
        function [fqfp,s,r] = t4imgs_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',          @this.lexist_4dfp);
            addRequired( ip, 'out',         @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.t4imgs_4dfp__( ...
                sprintf('%s %s %s', ...
                ip.Results.options, ...
                myfileprefix(ip.Results.in), ...
                myfileprefix(ip.Results.out)));
            fqfp = myfileprefix(ip.Results.out);
        end
        function      [s,r] = t4_ident(this, varargin)
            ip = inputParser;
            addRequired(ip, 't4file', @ischar);
            parse(ip, varargin{:});
            
            [s,r] = this.t4_ident__(ip.Results.t4file);
        end
        function   [t4,s,r] = t4_inv(this, varargin)
            ip = inputParser;
            addRequired( ip, 'in',          @(x) lexist(x, 'file'));
            addParameter(ip, 'out',         this.filenameT4inv(varargin{1}), ...
                                            @ischar);
            addParameter(ip, 'options', '', @ischar);
            parse(ip, varargin{:});
            t4 = ip.Results.out;
            
            [s,r] = this.t4_inv__(sprintf('%s %s %s', ip.Results.options, ip.Results.in, t4));
            try
                dbbash(sprintf('chmod 664 %s', ip.Results.out));
            catch ME
                handexcept(ME);
            end
        end
        function   [t4,s,r] = t4_mul(this, varargin)
            ip = inputParser;
            addRequired(ip, 'A2B', @(x) lexist(x, 'file'));
            addRequired(ip, 'B2C', @(x) lexist(x, 'file'));
            addOptional(ip, 'A2C', this.filenameT4mul(varargin{1:2}), ...
                                   @ischar);
            parse(ip, varargin{:});
            t4 = ip.Results.A2C;
            
            [s,r] = this.t4_mul__(sprintf('%s %s %s', ip.Results.A2B, ip.Results.B2C, t4));
        end
        function      [s,r] = t4_resolve(this, varargin)
            ip = inputParser;
            addRequired( ip, 'output',        @ischar);
            addRequired( ip, 'filenames',     @ischar);
            addParameter(ip, 'options',   '', @ischar);
            addParameter(ip, 'log',       '', @ischar);
            parse(ip, varargin{:});
            log = ip.Results.log;
            if (~isempty(ip.Results.log))
                log = [' &> ' log];
            end
            
            [s,r] = this.t4_resolve__( ...
                sprintf(' %s -o%s %s %s', ip.Results.options, ip.Results.output, ip.Results.filenames, log));
        end
        function      [s,r] = zero_slice_4dfp(this, varargin)
            ip = inputParser;
            addRequired( ip, 'fqfp',   @ischar);
            addRequired( ip, 'axis',   @ischar);
            addRequired(ip, 'index1',  @isnumeric);
            addRequired(ip, 'indexF',  @isnumeric);
            addOptional(ip, 'fdfpOut', sprintf('%s_%s%ito%i', varargin{1:4}), @ischar);
            parse(ip, varargin{:});            
            
            [s,r] = this.zero_slice_4dfp__( ...
                sprintf(' %s %s %i %i %s', ip.Results.fqfp, ip.Results.axis, ip.Results.index1, ip.Results.indexF, ip.Results.fdfpOut));
        end        
        
 		function this = FourdfpVisitor(varargin)
 			%% FOURDFPVISITOR
 			%  Usage:  this = FourdfpVisitor()

            if (this.ASSERT_PLATFORM)
                this.assertPlatform;
            end
        end        
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function [t4,fqfp,s,r] = ...
                       align_mode(this, varargin)
            %% ALIGN_MODE calls imgreg_4dfp with the specified mode and writes a log.
            %  @required mode is numeric; cf. ${TRANSFER}/imgreg_4dfp.txt.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   "; default := 'none'
            %  @param sourceMask "; default := 'none'
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param t4img_4dfp is logical.
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            ip = inputParser;
            addRequired( ip, 'mode',               @isnumeric);
            addParameter(ip, 'dest',       '',     @this.lexist_4dfp);
            addParameter(ip, 'source',     '',     @this.lexist_4dfp);
            addParameter(ip, 'destMask',   'none', @ischar);
            addParameter(ip, 'sourceMask', 'none', @ischar);
            addParameter(ip, 'destBlur',   0,      @isnumeric);
            addParameter(ip, 'sourceBlur', 0,      @isnumeric);
            addParameter(ip, 't40',        this.transverse_t4, @(x) lexist(x, 'file')); 
            addParameter(ip, 't4',         '',     @ischar);
            addParameter(ip, 't4img_4dfp', true,   @islogical);
            addParameter(ip, 'out',        '',     @ischar);
            addParameter(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});   
            dest       = ip.Results.dest;
            source     = ip.Results.source;
            destMask   = ip.Results.destMask;
            sourceMask = ip.Results.sourceMask;
            t4         = ip.Results.t4;
            log        = ip.Results.log;
            
            if (strcmp(dest, source))
                t4 = this.filenameT4(source, dest);
                copyfile(this.transverse_t4, t4); % identity transformation
                fqfp = '';
                return
            end
            if (isempty(t4))
                t4 = this.filenameT4(source, dest); end
            if (ip.Results.destBlur > 0)
                dest = this.imgblur_4dfp(dest, ip.Results.destBlur); end
            if (ip.Results.sourceBlur > 0)
                source = this.imgblur_4dfp(source, ip.Results.sourceBlur); end
            
            if (~strcmp(ip.Results.t40, t4))
                [s,r] = copyfile(ip.Results.t40, t4, 'f'); %#ok<ASGLU>
            end
            [s,r] = dbbash(sprintf('chmod 777 %s', t4)); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, ip.Results.mode, log);
            if (ip.Results.t4img_4dfp)
                [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'options', ['-O' ip.Results.dest]);
            else                
                [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'out', ip.Results.out, 'options', ['-O' ip.Results.dest]);
            end
        end
        function [t4,fqfp,s,r] = ...
                       align_modes2(this, varargin)
            %% ALIGN_MODES2 calls imgreg_4dfp twice successively with the specified mode and writes a log.
            %  @required mode is numeric; cf. ${TRANSFER}/imgreg_4dfp.txt.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   "; default := 'none'
            %  @param sourceMask "; default := 'none'
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param t4img_4dfp is logical.
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            ip = inputParser;
            addRequired( ip, 'mode',               @isnumeric);
            addParameter(ip, 'dest',       '',     @this.lexist_4dfp);
            addParameter(ip, 'source',     '',     @this.lexist_4dfp);
            addParameter(ip, 'destMask',   'none', @ischar);
            addParameter(ip, 'sourceMask', 'none', @ischar);
            addParameter(ip, 'destBlur',   0,      @isnumeric);
            addParameter(ip, 'sourceBlur', 0,      @isnumeric);
            addParameter(ip, 't40',        this.transverse_t4, @(x) lexist(x, 'file')); 
            addParameter(ip, 't4',         '',     @ischar);
            addParameter(ip, 't4img_4dfp', true, @islogical);
            addParameter(ip, 'out',        '',     @ischar);
            addParameter(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});   
            dest       = ip.Results.dest;
            source     = ip.Results.source;
            destMask   = ip.Results.destMask;
            sourceMask = ip.Results.sourceMask;
            t4         = ip.Results.t4;
            log        = ip.Results.log;
            
            if (strcmp(dest, source))
                t4 = this.filenameT4(source, dest);
                copyfile(this.transverse_t4, t4); % identity transformation
                fqfp = '';
                return
            end
            if (isempty(t4))
                t4 = this.filenameT4(source, dest); end
            if (ip.Results.destBlur > 0)
                dest = this.imgblur_4dfp(dest, ip.Results.destBlur); end
            if (ip.Results.sourceBlur > 0)
                source = this.imgblur_4dfp(source, ip.Results.sourceBlur); end
            
            if (~strcmp(ip.Results.t40, t4))
                [s,r] = copyfile(ip.Results.t40, t4, 'f'); %#ok<ASGLU>
            end
            [s,r] = dbbash(sprintf('chmod 777 %s', t4)); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, ip.Results.mode, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, ip.Results.mode, log); 
            if (ip.Results.t4img_4dfp)
                if (isempty(ip.Results.out))
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'options', ['-O' ip.Results.dest]);
                else
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'out', ip.Results.out, 'options', ['-O' ip.Results.dest]);
                end
            else                
                fqfp = ''; 
            end
        end
        function [t4,fqfp,s,r] = ...
                       align_modes5(this, varargin)
            %% ALIGN_MODES5 calls imgreg_4dfp with cross modes := [4099 4099 3075 2051 10243] 
            %  or common modes := [4099 4099 3075 2051 10243] + 256 and then writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   " ; default := 'none'
            %  @param sourceMask " ; default := 'none'
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest;   default := 0.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source; default := 0.
            %  @param t40        is the initial t4-file for the transformation; default := transverse.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param useCommonModal:     is logical, default := false; cf. ${TRANSFER}/cross_modal_*.ps.
            %  @param useMetricGradient:  is logical, default := true;  cf. ${TRANSFER}/imgreg_4dfp.txt.
            %  @param t4img_4dfp is logical.
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            ip = inputParser;
            addParameter(ip, 'dest',       '',     @this.lexist_4dfp);
            addParameter(ip, 'source',     '',     @this.lexist_4dfp);
            addParameter(ip, 'destMask',   'none', @ischar);
            addParameter(ip, 'sourceMask', 'none', @ischar);
            addParameter(ip, 'destBlur',   0,      @isnumeric);
            addParameter(ip, 'sourceBlur', 0,      @isnumeric);
            addParameter(ip, 't40',        this.transverse_t4, @(x) lexist(x, 'file')); 
            addParameter(ip, 't4',         '',     @ischar);
            addParameter(ip, 'useCommonModal',     false, @islogical);
            addParameter(ip, 'useMetricGradient',  true, @islogical);
            addParameter(ip, 't4img_4dfp',         true, @islogical);
            addParameter(ip, 'out',        '',     @ischar);
            addParameter(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});
            dest       = ip.Results.dest;
            source     = ip.Results.source;
            destMask   = ip.Results.destMask;
            sourceMask = ip.Results.sourceMask;
            t4         = ip.Results.t4;
            log        = ip.Results.log;
            if (ip.Results.useCommonModal)
                madj = 256;
                if (verbose)
                    fprintf('mlfourdfp.FourdfpVisitor.align_modes6:  using common modal methods for early steps'); end
            else
                madj = 0;
                if (verbose)
                    fprintf('mlfourdfp.FourdfpVisitor.align_modes6:  using cross modal registration'); end
            end
            
            if (strcmp(dest, source))
                t4 = this.filenameT4(source, dest);
                copyfile(this.transverse_t4, t4); % identity transformation
                fqfp = '';
                return
            end
            if (isempty(t4))
                t4 = this.filenameT4(source, dest); end
            if (ip.Results.destBlur > 0)
                dest = this.imgblur_4dfp(dest, ip.Results.destBlur); end
            if (ip.Results.sourceBlur > 0)
                source = this.imgblur_4dfp(source, ip.Results.sourceBlur); end
            
            if (~strcmp(ip.Results.t40, t4))
                [s,r] = copyfile(ip.Results.t40, t4, 'f'); %#ok<ASGLU>
            end
            [s,r] = dbbash(sprintf('chmod 777 %s', t4)); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, 'none',   source, 'none',     t4, 4099+madj, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, 'none',   source, 'none',     t4, 4099+madj, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 3075+madj, log); %#ok<ASGLU>   
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051+madj, log); %#ok<ASGLU>
            if (~ip.Results.useMetricGradient)
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051,      log); 
            else
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 10243,      log); 
            end
            if (ip.Results.t4img_4dfp)
                if (isempty(ip.Results.out))
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'options', ['-O' ip.Results.dest]);
                else
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'out', ip.Results.out, 'options', ['-O' ip.Results.dest]);
                end
            else
                fqfp = ''; 
            end
        end
        function [t4,fqfp,s,r] = ...
                       align_modes6(this, varargin)
            %% ALIGN_MODES6 calls imgreg_4dfp with cross modes := [4099 4099 3075 2051 2051 10243] 
            %  or common modes := [4099 4099 3075 2051 2051 10243] + 256 and then writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   " ; default := 'none'
            %  @param sourceMask " ; default := 'none'
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest;   default := 0.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source; default := 0.
            %  @param t40        is the initial t4-file for the transformation; default := transverse.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param useCommonModal:     is logical, default := false; cf. ${TRANSFER}/cross_modal_*.ps.
            %  @param useMetricGradient:  is logical, default := true;  cf. ${TRANSFER}/imgreg_4dfp.txt.
            %  @param t4img_4dfp is logical.
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            ip = inputParser;
            addParameter(ip, 'dest',       '',     @this.lexist_4dfp);
            addParameter(ip, 'source',     '',     @this.lexist_4dfp);
            addParameter(ip, 'destMask',   'none', @ischar);
            addParameter(ip, 'sourceMask', 'none', @ischar);
            addParameter(ip, 'destBlur',   0,      @isnumeric);
            addParameter(ip, 'sourceBlur', 0,      @isnumeric);
            addParameter(ip, 't40',        this.transverse_t4, @(x) lexist(x, 'file')); 
            addParameter(ip, 't4',         '',     @ischar);
            addParameter(ip, 'useCommonModal',     false, @islogical);
            addParameter(ip, 'useMetricGradient',  true, @islogical);
            addParameter(ip, 't4img_4dfp',         true, @islogical);
            addParameter(ip, 'out',        '',     @ischar);
            addParameter(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});
            dest       = ip.Results.dest;
            source     = ip.Results.source;
            destMask   = ip.Results.destMask;
            sourceMask = ip.Results.sourceMask;
            t4         = ip.Results.t4;
            log        = ip.Results.log;
            if (ip.Results.useCommonModal)
                madj = 256;
                if (verbose)
                    fprintf('mlfourdfp.FourdfpVisitor.align_modes6:  using common modal methods for early steps'); end
            else
                madj = 0;
                if (verbose)
                    fprintf('mlfourdfp.FourdfpVisitor.align_modes6:  using cross modal registration'); end
            end
            
            if (strcmp(dest, source))
                t4 = this.filenameT4(source, dest);
                copyfile(this.transverse_t4, t4); % identity transformation
                fqfp = '';
                return
            end
            if (isempty(t4))
                t4 = this.filenameT4(source, dest); end
            if (ip.Results.destBlur > 0)
                dest = this.imgblur_4dfp(dest, ip.Results.destBlur); end
            if (ip.Results.sourceBlur > 0)
                source = this.imgblur_4dfp(source, ip.Results.sourceBlur); end
            
            if (~strcmp(ip.Results.t40, t4))
                [s,r] = copyfile(ip.Results.t40, t4, 'f'); %#ok<ASGLU>
            end
            [s,r] = dbbash(sprintf('chmod 777 %s', t4)); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, 'none',   source, 'none',     t4, 4099+madj, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, 'none',   source, 'none',     t4, 4099+madj, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 3075+madj, log); %#ok<ASGLU>   
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051+madj, log); %#ok<ASGLU>
            if (~ip.Results.useMetricGradient)
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051,      log); 
            else
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 10243+madj, log); %#ok<ASGLU>
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 10243,      log); 
            end
            if (ip.Results.t4img_4dfp)
                if (isempty(ip.Results.out))
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'options', ['-O' ip.Results.dest]);
                else
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'out', ip.Results.out, 'options', ['-O' ip.Results.dest]);
                end
            else
                fqfp = ''; 
            end
        end
        function [t4,fqfp,s,r] = ...
                       align_modes7(this, varargin)
            %% ALIGN_MODES7 calls imgreg_4dfp with cross modes := [4099 4099 3075 3075 2051 2051 10243] 
            %  or common modes := [4099 4099 3075 3075 2051 2051 10243] + 256 and then writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   " ; default := 'none'
            %  @param sourceMask " ; default := 'none'
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest;   default := 0.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source; default := 0.
            %  @param t40        is the initial t4-file for the transformation; default := transverse.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param useCommonModal:     is logical, default := false; cf. ${TRANSFER}/cross_modal_*.ps.
            %  @param useMetricGradient:  is logical, default := true;  cf. ${TRANSFER}/imgreg_4dfp.txt.
            %  @param t4img_4dfp is logical.
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            ip = inputParser;
            addParameter(ip, 'dest',       '',     @this.lexist_4dfp);
            addParameter(ip, 'source',     '',     @this.lexist_4dfp);
            addParameter(ip, 'destMask',   'none', @ischar);
            addParameter(ip, 'sourceMask', 'none', @ischar);
            addParameter(ip, 'destBlur',   0,      @isnumeric);
            addParameter(ip, 'sourceBlur', 0,      @isnumeric);
            addParameter(ip, 't40',        this.transverse_t4, @(x) lexist(x, 'file')); 
            addParameter(ip, 't4',         '',     @ischar);
            addParameter(ip, 'useCommonModal',     false, @islogical);
            addParameter(ip, 'useMetricGradient',  true, @islogical);
            addParameter(ip, 't4img_4dfp',         true, @islogical);
            addParameter(ip, 'out',        '',     @ischar);
            addParameter(ip, 'log',        '/dev/null', @ischar);
            parse(ip, varargin{:});
            dest       = ip.Results.dest;
            source     = ip.Results.source;
            destMask   = ip.Results.destMask;
            sourceMask = ip.Results.sourceMask;
            t4         = ip.Results.t4;
            log        = ip.Results.log;
            if (ip.Results.useCommonModal)
                madj = 256;
                if (verbose)
                    fprintf('mlfourdfp.FourdfpVisitor.align_modes6:  using common modal methods for early steps'); end
            else
                madj = 0;
                if (verbose)
                    fprintf('mlfourdfp.FourdfpVisitor.align_modes6:  using cross modal registration'); end
            end
            
            if (strcmp(dest, source))
                t4 = this.filenameT4(source, dest);
                copyfile(this.transverse_t4, t4); % identity transformation
                fqfp = '';
                return
            end
            if (isempty(t4))
                t4 = this.filenameT4(source, dest); end
            if (ip.Results.destBlur > 0)
                dest = this.imgblur_4dfp(dest, ip.Results.destBlur); end
            if (ip.Results.sourceBlur > 0)
                source = this.imgblur_4dfp(source, ip.Results.sourceBlur); end
            
            if (~strcmp(ip.Results.t40, t4))
                [s,r] = copyfile(ip.Results.t40, t4, 'f'); %#ok<ASGLU>
            end
            [s,r] = dbbash(sprintf('chmod 777 %s', t4)); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, 'none',   source, 'none',     t4, 4099+madj, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, 'none',   source, 'none',     t4, 4099+madj, log); %#ok<ASGLU>
            [s,r] = this.imgreg_4dfp(dest, destMask, source, 'none',     t4, 3075+madj, log); %#ok<ASGLU>  
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 3075+madj, log); %#ok<ASGLU>   
            [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051+madj, log); %#ok<ASGLU>
            if (~ip.Results.useMetricGradient)
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051+madj, log); %#ok<ASGLU>
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 2051,      log); 
            else
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 10243+madj, log); %#ok<ASGLU>
                [s,r] = this.imgreg_4dfp(dest, destMask, source, sourceMask, t4, 10243,      log); 
            end
            if (ip.Results.t4img_4dfp)
                if (isempty(ip.Results.out))
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'options', ['-O' ip.Results.dest]);
                else
                    [fqfp,s,r] = this.t4img_4dfp(t4, ip.Results.source, 'out', ip.Results.out, 'options', ['-O' ip.Results.dest]);
                end
            else
                fqfp = ''; 
            end            
        end
        function       assertPlatform(this)
            if (lgetenv('DEBUG'))
                return
            end            
            assert(strcmp(computer, 'GLNXA64'));
            [~,hn] = mlbash('hostname -f');
            assert(lstrfind(hn, this.FOURDFP_HOSTS));
        end        
        function img = normalizeFramesBySums(~, img)
            assert(length(size(img)) == 4);
            for w = 1:size(img, 4)
                img(:,:,:,w) = img(:,:,:,w) / dipmean(img(:,:,:,w));
            end
        end
    end

    %% PRIVATE
    
    methods (Static, Access = private)
        function s = safesprintf(unsafeStr, fp, idxs)
            idx1 = idxs(1);
            idx2 = idx1;
            try
                switch (unsafeStr)
                    case '_on_'
                        idx2 = idx1+4;
                        s = sprintf('%sOn%s%s', fp(1:idx1-1), upper(fp(idx2)), fp(idx2+1:end));
                    case '_op_'
                        idx2 = idx1+4;
                        s = sprintf('%sOp%s%s', fp(1:idx1-1), upper(fp(idx2)), fp(idx2+1:end));
                    case '+'
                        idx2 = idx1+1;
                        s = sprintf('%s%s%s', fp(1:idx1-1), upper(fp(idx2)), fp(idx2+1:end));
                    case '.a2009s'
                        idx2 = idx1+1;
                        s = sprintf('%s%s%s', fp(1:idx1-1), upper(fp(idx2)), fp(idx2+1:end));
                    otherwise
                        error('mlfourdfp:unsupportedSwitchcase', 'FourdfpVisitor.safesprintf.fp->%s', fp)
                end
            catch ME
                assert(idx1 > 1, ...
                    'FourdfpVisitor.safesprintf.fp->%s may be missing prefix', fp);
                assert(length(fp) >= idx2, ...
                    'FourdfpVisitor.safesprintf.fp->%s may be missing suffix', fp);
                handexcept(ME);
            end
        end
    end
    
    methods (Access = private)
        function [s,r] = CT2mpr_4dfp__(~, args)
            %% CT2MPR_4DFP
            % $Id: CT2mpr_4dfp,v 1.3 2016/02/11 22:45:40 avi Exp $
            % Usage:	CT2mpr_4dfp <4dfp mprage> <4dfp CT> [options]
            % e.g.:	CT2mpr_4dfp PT34_mpr1 ../CT090202/OTSN80_CT -T/data/cninds01/data2/atlas/TRIO_Y_NDC
            % 	options
            % 	-T<target>	specify atlas target (<target> may include absolute path)
            % 	-m	run additional final registration using a fat MP-RAGE mask
            % 	setecho	set echo
            % 	debug	debug mode
            % 	redo	recompute (e.g., after manual t4file adjustment)
            % N.B.:	CT2mpr_4dfp assumes that <4dfp mprage> is in the current working directory
            % 	and that its atlas transform, e.g., PT34_mpr1_to_TRIO_Y_NDC_t4
            % 	exists and is in the current working directory     
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('CT2mpr_4dfp %s', args));       
        end
        function [s,r] = C2T_4dfp__(~, args)
            %% C2T_4dfp__
            % $Id: C2T_4dfp.c,v 1.6 2007/08/31 05:03:24 avi Exp $
            % Usage:	C2T_4dfp <(4dfp) imgroot> [(4dfp) outroot]
            %  e.g.,	C2T_4dfp vm6c_mpr
            %  e.g.,	C2T_4dfp vm6c_mpr vm6c_mprT
            % 	option
            % 	-@<b|l>	output big or little endian (default input endian)
            % N.B.:	default output root = <imgroot>"T"
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('C2T_4dfp %s', args));
        end
        function [s,r] = S2T_4dfp__(~, args)
            %% S2T_4dfp__
            % $Id: S2T_4dfp.c,v 1.6 2007/08/31 05:03:24 avi Exp $
            % Usage:	S2T_4dfp <(4dfp) imgroot> [(4dfp) outroot]
            %  e.g.,	S2T_4dfp vm6c_mpr
            %  e.g.,	S2T_4dfp vm6c_mpr vm6c_mprT
            % 	option
            % 	-@<b|l>	output big or little endian (default input endian)
            % N.B.:	default output root = <imgroot>"T"
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('S2T_4dfp %s', args));
        end
        function [s,r] = T2C_4dfp__(~, args)
            %% T2C_4dfp__
            % $Id: T2C_4dfp.c,v 1.6 2007/08/31 05:03:24 avi Exp $
            % Usage:	T2C_4dfp <(4dfp) imgroot> [(4dfp) outroot]
            %  e.g.,	T2C_4dfp vm6c_mpr
            %  e.g.,	T2C_4dfp vm6c_mpr vm6c_mprT
            % 	option
            % 	-@<b|l>	output big or little endian (default input endian)
            % N.B.:	default output root = <imgroot>"T"
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('T2C_4dfp %s', args));
        end
        function [s,r] = T2S_4dfp__(~, args)
            %% T2S_4dfp__
            % $Id: S2T_4dfp.c,v 1.6 2007/08/31 05:03:24 avi Exp $
            % Usage:	T2S_4dfp <(4dfp) imgroot> [(4dfp) outroot]
            %  e.g.,	T2S_4dfp vm6c_mpr
            %  e.g.,	T2S_4dfp vm6c_mpr vm6c_mprT
            % 	option
            % 	-@<b|l>	output big or little endian (default input endian)
            % N.B.:	default output root = <imgroot>"T"
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('T2S_4dfp %s', args));
        end
        function [s,r] = IFhdr_to_4dfp__(~, args)
            %% IFHDR_TO_4DFP__
            % $Id: IFhdr_to_4dfp,v 1.1 2014/02/25 suy Exp $
            % Usage: IFhdr_to_4dfp IFfstr outroot

            assert(ischar(args));
            [s,r] = dbbash(sprintf('IFhdr_to_4dfp %s', args));
        end
        function [s,r] = actmapf_4dfp__(~, args)
            %% ACTMAPF_4DFP__
            % $Id: actmapf_4dfp.c,v 1.34 2016/02/12 00:59:43 avi Exp $
            % Usage:	actmapf_4dfp <format|fmtfile> <4dfp|conc input>
            %  e.g.,	actmapf_4dfp -zu "3x3(11+4x15-)" b1_rmsp_dbnd_xr3d_norm
            %  e.g.,	actmapf_4dfp -aanatomy -c10 -u "+" ball_dbnd_xr3d.conc
            %  e.g.,	actmapf_4dfp -zu "4x124+" b1_rmsp_dbnd_xr3d -wweights.txt
            % 	        option
            % 	        -a<str>	specify 4dfp output root trailer (default = "actmap")
            % 	        -c<flt>	scale output by specified factor
            % 	        -u	scale weights to unit variance
            % 	        -z	adjust weights to zero sum
            % 	        -R	compute relative modulation (default absolute)
            % 	        -w<weight file>	read (text) weights from specified filename
            % 	        -@<b|l>	output big or little endian (default input endian)
            % N.B.:	    conc files must have extension "conc"
            % N.B.:	    when using weight files 'x' frames in format are not counted
            % N.B.:	    relative modulation images are zeroed where mean intensity < 0.5*whole_image_mode
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('actmapf_4dfp %s', args));
        end
        function [s,r] = compute_defined_4dfp__(~, args)
            %% COMPUTE_DEFINED_4DFP__
            % $Id: compute_defined_4dfp.c,v 1.8 2016/02/14 06:37:36 avi Exp $
            % Usage:	compute_defined_4dfp <4dfp|conc input>
            % 	option
            % 	-z	count zero voxels as undefined (default defined)
            % 	-f<str>	specify frames-to-count format (default count all)
            % 	-F<str>	read frames-to-count format from specified file

            assert(ischar(args));
            [s,r] = dbbash(sprintf('compute_defined_4dfp %s', args));
        end
        function [s,r] = copy_4dfp__(~, args)
            %% COPY_4DFP__
            % Usage:	copy_4dfp source destination
            %  e.g.,	copy_4dfp bad_4dfp_name better_4dfp_name
            %  
            % copy_4dfp works only with 4dfp files; it does not work on directories
            assert(ischar(args));
            [s,r] = dbbash(sprintf('copy_4dfp %s', args));
        end
        function [s,r] = cropfrac_4dfp__(~, args)
            %% CROPFRAC_4DFP__
            % Usage:/mnt/hgfs/Local/bin/cropfrac_4dfp <fraction to crop> <4dfp to crop> <4dfp cropped>
            % e.g.,/mnt/hgfs/Local/bin/cropfrac_4dfp 0.5 NP995_09ho1_v1 ho1v1
            % N.B.:crops only along x,y axes; cf. crop_4dfp for details

            assert(ischar(args));
            [s,r] = dbbash(sprintf('cropfrac_4dfp %s', args));
        end
        function [s,r] = crop_4dfp__(~, args)
            %% CROP_4DFP__
            % $Id: crop_4dfp.c,v 1.10 2010/11/12 05:54:27 avi Exp $
            % Usage:	crop_4dfp <(4dfp) inroot> [(4dfp) outroot]
            % 	        option
            % 	        -<x|y|z><int>[to[<int>]	specify x y z crop limits
            % 	        -s<x|y|z><int>	scroll specified axis by specified number of pixels
            % 	        -f	interpret specifications under 4dfp<->analyze flips
            % 	        -Z	zero voxels instead of physically cropping
            % 	        -@<b|l>	output big or little endian (default input endian)
            % N.B.:	    crop limit indices count from 1
            % N.B.:	    scrolling is done after cropping
            % N.B.:	    default (4dfp) output root is <(4dfp) inroot>"_crop"
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('crop_4dfp %s', args));
        end
        function [s,r] = dcm_sort__(~, args)
            %% DCM_SORT__
            % $Id: dcm_sort,v 1.21 2015/07/24 06:21:38 avi Exp $
            % Usage:	dcm_sort <dicom_directory>
            % e.g.,	dcm_sort /data/petsun52/data1/JHILL/04271737
            % e.g.,	dcm_sort /cdrom/botv/10251349 -p930589002 -c
            % 	options
            % 	-d	verbose debug mode
            % 	-c	copy files (default symbolically link)
            % 	-t	toggle use of -t in call to dcm_dump_file (default ON)
            % 	-i	take files with integer filenames
            % 	-e<ext>	take files with specified extension
            % 	-r<str>	take files with filenames containing specified string (default is 7 digits)
            % 	-p<str>	take files only with dicom field 'PAT Patient Name' matching specified string
            % N.B.:	dcm_sort removes existing single study subdirectories
            % N.B.:	dcm_sort puts unclassifiable DICOMs into subdirectory study0

            assert(ischar(args));
            [s,r] = dbbash(sprintf('dcm_sort %s', args));
        end
        function [s,r] = dcm_sort_PPG__(~, args)
            %% DCM_SORT_PPG__

            assert(ischar(args));
            [s,r] = dbbash(sprintf('dcm_sort_PPG.csh %s', args));
        end
        function [s,r] = dcm_to_4dfp__(~, args)
            %% DELETE_4DFP__
            % $Id: 507d5bd9b20fe72fe5ccd013aad04875baf65952 $
            % Usage:	dcm_to_4dfp [-b base] [-d gggg eeee] [-f] [-g] [-u] file(s)
            % Slice Spacing Options: [-c] [-t <flt> or S or T]
            % Slice Position Options: [-X] [-Y] [-Z]
            %  e.g.,	dcm_to_4dfp *
            %    or,	dcm_to_4dfp -b ID101 -f -g -u *IMA
            %    or,	dcm_to_4dfp -d 0008 0030 -t 4.98 -g *.dcm
            %    or,	dcm_to_4dfp -b P0089 -t T -g mydir/*
            % Options:
            % 	[-b base] Output base filename follows the -b.
            % 	[-c]	    Slice Spacing: By Image Position (0020 0032).
            % 	[-d gggg eeee] Divide series by group and element number.
            % 		 ** Default will divide volumes using ID series time (0008 0031).
            % 	[-f]	    Directories will be created, and dicom files will be moved.
            % 	[-g]	    Add image name, XYZ relative position, and number to rec file.
            %
            % 	[-q]       Slice Spacing: Do not compute by Image Position.
            % 	[-r]       Rescale: Use the rescale slope and intercept fields.
            % 	[-t <flt>] Slice Spacing: Use input value.[-t <flt>]
            % 	[-t T]     Slice Spacing: Use Slice Thickness 0018 0050.[-t T]
            % 	[-t S]     Slice Spacing: Use Slice Spacing 0018 0088 [-t S](** Default)
            % 	[-u]	Output files named using sequence tag 0018 0024 plus number.
            %
            % 		4dfp Coordinant System is determined by Image Position (0020 0032).
            % 		Multivolume and BOLD images are ordered by REL Image Number (0020 0013).
            % 	[-X]	Sagittal:	image positions will be ordered low to high
            % 	[-Y]	Coronal:	image positions will be high to low
            % 	[-Z]	Transverse:	image positions will be high to low
            % 		** Default is transverse ordered by REL Image Number (0020 0013).
            % 	[-@ <b|l>]	output big or little endian (default CPU endian)

            assert(ischar(args));
            [s,r] = dbbash(sprintf('dcm_to_4dfp %s', args));
        end
        function [s,r] = delete_4dfp__(~, args)
            %% DELETE_4DFP__
            %  Usage:	delete_4dfp 4dfp-basename
            %  e.g.,	delete_4dfp files_to_delete # .4dfp.*
 
            assert(ischar(args));
            [s,r] = dbbash(sprintf('delete_4dfp %s', args));
        end
        function bn    = ensureImgregFiles(this, fqfp)
            if (strcmp('none', fqfp))
                bn = fqfp;
                return
            end
            bn = mybasename(fqfp);
            s = 0; r = '';
            if (~lexist([bn '.4dfp.img']))
                try
                    [s,r] = this.lns_4dfp(fqfp);  %#ok<ASGLU>
                catch ME
                    fprintf('s->%i, r->%s\n', s, r); 
                    handerror(ME);
                end
            end
        end
        function [s,r] = extract_frame_4dfp__(~, args)
            %% EXTRACT_FRAME_4DFP__
            % usage:	extract_frame_4dfp <(4dfp) stack> <(int) frame>
            % e.g.:	    extract_frame_4dfp CDR.5to1+ 3
            % 	        option
            % 	        -o<str>	specifiy output 4dfp fileroot (default = <stack>_frame<frame>)
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('extract_frame_4dfp %s', args));
        end
        function [s,r] = freesurfer2mpr_4dfp__(~, args)
            %% FREESURFER2MPR_4DFP__
            % $Id: freesurfer2mpr_4dfp,v 1.9 2017/08/17 22:27:00 avi Exp $
            % Usage:  freesurfer2mpr_4dfp <(4dfp) mpr> <(4dfp) orig> [options]
            % e.g.,   freesurfer2mpr_4dfp vc1234_654-3[.4dfp.img] vc1234_orig
            % e.g.,   freesurfer2mpr_4dfp vc1234_654-3 vc1234_orig -T711-2V apply
            %         options
            %         -skew           general affine orig->mpr registeration (default 6 parameter rigid body)
            %         -T<target>      specify atlas representative target
            %         -a<segimg>      add named (4dfp format) freesurfer segemntation result to "apply" list
            %         apply   proceed directly to transform (4dfp format) segmentations
            %         force   force atlas transformation of segmentation results even if it already exists
            %         setecho set echo
            % N.B.:   <(4dfp) orig> is the freesurfer-resampled 256x256x256 coronal mpr
            % N.B.:   the default "apply" list includes (4dfp format) images named *parc* and *aseg*
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('freesurfer2mpr_4dfp %s', args));
        end
        function [s,r] = flip_4dfp__(~, args)
            %% FLIP_4DFP__
            % $Id: flip_4dfp.c,v 1.6 2007/05/04 00:36:21 avi Exp $
            % Usage:	flip_4dfp <(4dfp) image> [(4dfp) output]
            % e.g.,	flip_4dfp -yz vc345 vc345_flipyz
            % 	option
            % 	-x	flip x
            % 	-y	flip y
            % 	-z	flip z
            % 	-@<b|l>	output big or little endian (default input endian)
            % N.B.:	default output fileroot = <image>_flip[xyz]
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('flip_4dfp %s', args));
        end
        function [s,r] = gauss_4dfp__(~, args)
            %% GAUSS_4DFP__
            % $Id: gauss_4dfp.c,v 1.16 2007/05/04 01:26:56 avi Exp $
            % Usage:	gauss_4dfp <4dfp|conc input> f_half [outroot]
            %  e.g.,	gauss_4dfp pt349_study9to9 0.1
            %  e.g.,	gauss_4dfp p1234ho5 0.7 p1234ho5_g7
            % 	        options
            % 	        -@<b|l>	output big or little endian (default input endian)
            % 	        -w	(wrap) suppress x and y padding
            % 	        -d	differentiate
            % N.B.:	    f_half is half frequency in 1/cm
            % N.B.:	    default output root is <inroot>_g<10*f_half>
            % N.B.:	    FWHM*f_half = (2ln2/pi) = 0.4412712
            % N.B.:	    conc files must have extension "conc"
            % N.B.:	    user outroot specification not possible with conc files
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('gauss_4dfp %s', args));
        end
        function [s,r] = imgblur_4dfp__(~, args)
            %% IMGBLUR_4DFP__
            % $Id: imgblur_4dfp.c,v 1.9 2009/07/18 01:30:43 avi Exp $
            % Usage:	imgblur_4dfp [options] <image_file> <FWHM_in_mm>
            %  e.g.,	imgblur_4dfp -yz vc345 5.5
            % 	        option
            % 	        -x	slective x blur
            % 	        -y	slective y blur
            % 	        -z	slective z blur
            % 	        -@<b|l>	output big or little endian (default input endian)
            % N.B.:	    default blur is 3D isotropic
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('imgblur_4dfp %s', args));
        end
        function [s,r] = imgopr_4dfp__(~, args)
            %% IMGOPR_4DFP__
            % $Id: imgopr_4dfp.c,v 1.14 2010/12/23 02:43:21 avi Exp $
            % Usage:	imgopr_4dfp -<operation><(4dfp) outroot> <(4dfp) image1> <(4dfp) image2> ...
            % 	operation
            % 	a	add
            % 	s	subtract (image1 - image2)
            % 	p	product
            % 	r	ratio (image1 / image2)
            % 	e	mean (expectation)
            % 	v	variance
            % 	g	geometric mean
            % 	n	count defined (see -u option) voxels
            % 	x	voxelwize maximum
            % 	y	voxelwize minimum
            % 	G	report serial number (counting from 1) of image with greatest value
            % 	P	unsplit multiple ROIs into fidl compatible ROI file
            % 	option
            % 	-u	count only defined (not NaN or 1.e-37 or 0.0) voxels
            % 	-R	suppress creation of rec file
            % 	-N	output undefined voxels as NaN
            % 	-Z	output undefined voxels as 0
            % 	-E	output undefined voxels as 1.E-37 (default)
            % 	-c<flt>	multiply output by specified scaling factor
            % 	-l<lst>	read input file names from specified list file
            % 	-@<b|l>	output big or little endian (default first input endian)
            % N.B.:	image dimensions must match except for binary operations {aspr} in which
            % 	a 1 volume second image may be paired with a multi-volume first image
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('imgopr_4dfp %s', args));
        end
        function [s,r] = imgreg_4dfp__(~, args)
            %% IMGREG_4DFP__
            % $Id: imgreg_4dfp.c,v 1.15 2007/08/10 03:27:30 avi Exp $
            % Usage:	imgreg_4dfp target_imag target_mask source_imag source_mask t4file mode
            % or:	    imgreg_4dfp target_imag        none source_imag source_mask t4file mode
            % or:	    imgreg_4dfp target_imag        none source_imag        none t4file mode
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('imgreg_4dfp %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.imgreg_4dfp__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = maskimg_4dfp__(~, args)
            %% MASKIMG_4DFP
            % $Id: maskimg_4dfp.c,v 1.17 2012/06/27 22:33:19 avi Exp $
            % Usage:	maskimg_4dfp <(4dfp) imgfile> <(4dfp) mskfile> <(4dfp) outfile>
            %  e.g.:	maskimg_4dfp -t23.2 va1234_mpr mask va1234_mpr_msk
            % 	option
            % 	-N	replace NaN in <imgfile> with corresponding <mskfile> value
            % 	-e	report to stdout mean <imgfile> within-mask value
            % 	-1	apply first frame of <mskfile> to all frames of <imgfile>
            % 	-R	suppress creation of rec file
            % 	-v<flt>	specify <outfile> uniform within-mask value
            % 	-p<flt>	specify <mskfile> threshold as percent of <mskfile> max
            % 	-t<flt>	specify <mskfile> threshold directly (default = 0.0)
            % 	-A	threshold mask by absolute value of <mskfile>
            % 	-@<b|l>	output big or little endian (default <imgfile> endian)
            % N.B.:	<imgfile> and <mskfile> may be the same

            assert(ischar(args));
            [s,r] = dbbash(sprintf('maskimg_4dfp %s', args));
        end
        function [s,r] = move_4dfp__(~, args)
            %% MOVE_4DFP__
            % Usage: move_4dfp source destination
            % e.g.,  move_4dfp bad_4dfp_name better_4dfp_name
            %  
            % move_4dfp works only with 4dfp files; it does not work on directories
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('move_4dfp %s', args));
        end
        function [s,r] = mpr2atl_4dfp__(~, args)
            %% MPR2ATL_4DFP__
            % $Id: mpr2atl_4dfp,v 1.27 2014/02/16 02:46:17 avi Exp $
            % Usage:	mpr2atl_4dfp <mpr_anat> [options]
            % e.g.,	    mpr2atl_4dfp vc1234_654-3[.4dfp.img]
            % e.g.,	    mpr2atl_4dfp vc1234_654-3[.4dfp.img] -T/data/petsun23/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % 	        options
            % 	        711-2<C|O|Y|K|L|G|H|V|F>	specify 711-2? series atlas representative target image
            % 	        -T<target including path>	specify arbitrary     atlas representative target image
            % 	        -S<atlas space>			specify atlas space (default=711-2B space)
            % 	        crossmodal	use cross-modal mpr->target registration
            % 	        useold		suppress recomputation  of existing t4 file
            % 	        redo		suppress initialization of existing t4 file
            % 	        setecho		set echo
            % N.B.:	    <mpr_anat> may include a path, e.g., /data/petmr1/data7/stem9/scout/654-3
            % N.B.:	    <mpr_anat> must be in either ANALYZE short int or 4dfp format; ANALYZE will be converted to 4dfp
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('mpr2atl_4dfp %s', args));
        end
        function [s,r] = msktgen_4dfp__(~, args)
            %% MSKTGEN_4DFP__
            % $Id: msktgen_4dfp,v 1.21 2010/07/27 01:52:58 avi Exp $
            % Usage:	msktgen_4dfp <(4dfp) image> [threshold] -T<target including path>  -S<atlas space>
            % e.g.,	    msktgen_4dfp 4859-5_mpr
            % e.g.,	    msktgen_4dfp 4859-5_mpr -T/data/petsun29/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % N.B.:	    default threshold = 200
            % 	        Specifiy a higher threshold, e.g., 400 for a tighter mask and vice versa
            % N.B.:	    msktgen_4dfp uses the first legitimate atlas transform t4 file it sees in 
            % 	        the current working directory, i.e., one of <image>_to_711-2*_t4
            %  	        or  one of <image>_to_<target>_t4
            % N.B.:	-S specifies the atlas space to use. Atlas supported currently is 711-2B. -S must be used with -T option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('msktgen_4dfp %s', args));
        end
        function [s,r] = msktgen_b110_4dfp__(~, args)
            %% MSKTGEN_B110_4DFP__
            % $Id: msktgen_4dfp,v 1.21 2010/07/27 01:52:58 avi Exp $
            % Usage:	msktgen_4dfp <(4dfp) image> [threshold] -T<target including path>  -S<atlas space>
            % e.g.,	    msktgen_4dfp 4859-5_mpr
            % e.g.,	    msktgen_4dfp 4859-5_mpr -T/data/petsun29/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % N.B.:	    default threshold = 200
            % 	        Specifiy a higher threshold, e.g., 400 for a tighter mask and vice versa
            % N.B.:	    msktgen_4dfp uses the first legitimate atlas transform t4 file it sees in 
            % 	        the current working directory, i.e., one of <image>_to_711-2*_t4
            %  	        or  one of <image>_to_<target>_t4
            % N.B.:	-S specifies the atlas space to use. Atlas supported currently is 711-2B. -S must be used with -T option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('msktgen_b110_4dfp %s', args));
        end
        function [s,r] = msktgen2_4dfp__(~, args)
            %% MSKTGEN2_4DFP__
            % $Id: msktgen_4dfp,v 1.21 2010/07/27 01:52:58 avi Exp $
            % Usage:	msktgen_4dfp <(4dfp) image> [threshold] -T<target including path>  -S<atlas space>
            % e.g.,	    msktgen_4dfp 4859-5_mpr
            % e.g.,	    msktgen_4dfp 4859-5_mpr -T/data/petsun29/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % N.B.:	    default threshold = 200
            % 	        Specifiy a higher threshold, e.g., 400 for a tighter mask and vice versa
            % N.B.:	    msktgen_4dfp uses the first legitimate atlas transform t4 file it sees in 
            % 	        the current working directory, i.e., one of <image>_to_711-2*_t4
            %  	        or  one of <image>_to_<target>_t4
            % N.B.:	-S specifies the atlas space to use. Atlas supported currently is 711-2B. -S must be used with -T option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('msktgen2_4dfp %s', args));
        end
        function [s,r] = msktgen3_4dfp__(~, args)
            %% MSKTGEN3_4DFP__
            % $Id: msktgen_4dfp,v 1.21 2010/07/27 01:52:58 avi Exp $
            % Usage:	msktgen_4dfp <(4dfp) image> [threshold] -T<target including path>  -S<atlas space>
            % e.g.,	    msktgen_4dfp 4859-5_mpr
            % e.g.,	    msktgen_4dfp 4859-5_mpr -T/data/petsun29/data1/atlas/NP345_111[.4dfp.img] -S711-2B
            % N.B.:	    default threshold = 200
            % 	        Specifiy a higher threshold, e.g., 400 for a tighter mask and vice versa
            % N.B.:	    msktgen_4dfp uses the first legitimate atlas transform t4 file it sees in 
            % 	        the current working directory, i.e., one of <image>_to_711-2*_t4
            %  	        or  one of <image>_to_<target>_t4
            % N.B.:	-S specifies the atlas space to use. Atlas supported currently is 711-2B. -S must be used with -T option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('msktgen3_4dfp %s', args));
        end
        function [s,r] = nifti_4dfp__(~, args)
            % NIFTI_4DFP__
            % $Id: nifti_4dfp.c,v 1.9 2011/09/13 03:40:37 avi Exp $
            % Usage: nifti_4dfp <-4 or -n> <infile> <outfile> [options]
            %  e.g.: nifti_4dfp -n time_BOXzstat_333_t88.4dfp.ifh time_BOXzstat_333_t88.nii
            % 	options
            % 	-T <t4 file>	specify a t4 file to use converting TO NIfTI from 4dfp
            % 	-n	convert TO NIfTI from 4dfp
            % 	-4	convert TO 4dfp from NIfTI
            % 	-N	suppress saving of mmppix and center fields in output ifh
            % 	-@<val>	specify endianness for output, b or B for big, l or L for little
            % N.B.:	exactly one of -4 or -n must be specified
            % N.B.:	".4dfp.ifh" or ".nii" are appended to filenames specified without extension
            % N.B.:	option -N has effect only on converting nii->4dfp
            % N.B.:	option -T has effect only on converting 4dfp->nii
            assert(ischar(args));
            if (lstrcmp(computer, 'MACI64'))
                [s,r] = dbbash(sprintf('%s/Local/bin/nifti_4dfp %s', getenv('HOME'), args));
                return
            end
            [s,r] = dbbash(sprintf('nifti_4dfp %s', args));
        end
        function [s,r] = paste_4dfp__(~, args)
            %% PASTE_4DFP__
            % $Id: paste_4dfp.c,v 1.9 2008/10/17 02:53:48 avi Exp $
            % Usage:	paste_4dfp <inlist> <outfile>
            % 	        option
            % 	        -a	append successive epochs (default average)
            % 	        -p<int>	specify period in frames (default=1)
            % 	        -@<b|l>	output big or little endian (default initial input endian)
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('paste_4dfp %s', args));
        end
        function [s,r] = pseudo_dcm_sort__(~, args)
            %% PSEUDO_DCM_SORT__
            % usage:	pseudo_dcm_sort.csh <dicom directory>
            % e.g.:	pseudo_dcm_sort.csh RAW
            % N.B.:	dicom subdirectories must be numeric
            % 	option
            % 	-d	debug mode (set echo)
            % 	-s	DICOM files are within a subdirectory of numeric subdirectories (default directly in numeric subdirectory)
            % 	-e	identify DICOM files by specified extension (default extension = dcm)
            % 	-r	identify DICOM files by specified root (default root = 'MR*')
            % 	-i	take DICOM files with integer filenames
            % 	-t	toggle off use of -t in call to dcm_dump_file (default ON)
            % N.B.:	default subdirectory of numeric subdirectory is 'DICOM'

            assert(ischar(args));
            [s,r] = dbbash(sprintf('pseudo_dcm_sort.csh %s', args));
        end
        function [s,r] = scale_4dfp__(~, args)
            %% SCALE_4DFP__
            % $Id: scale_4dfp.c,v 1.14 2007/07/05 05:01:52 avi Exp $
            % Usage:	scale_4dfp <image_4dfp> <scale_factor> [options]
            % 	option
            % 	-E	preserve 1.0e-37 values (fidl NaN convention)
            % 	-a<str>	append trailer to output file name
            % 	-b<flt>	add specified constant to each voxel
            % 	-@<b|l>	output big or little endian (default input endian)
            % e.g.,	scale_4dfp b2_xfrm_avg 12
            % e.g.,	scale_4dfp b2_xfrm_avg 12 -b5 -ax12+5
            % N.B.:	<image_4dfp> is overwritten unless the trailer option is used
            % N.B.:	<scale_factor> must be specified for proper operation

            assert(ischar(args));
            [s,r] = dbbash(sprintf('scale_4dfp %s', args));
        end
        function [s,r] = sqrt_4dfp__(~, args)
            %% SQRT_4DFP__
            % $Id: sqrt_4dfp.c,v 1.7 2008/03/14 02:24:00 avi Exp $
            % Usage:	sqrt_4dfp <(4dfp) image> [outroot]
            % e.g.,	sqrt_4dfp vce20_mpr
            % 	-@<b|l>	output big or little endian (default input endian)
            % 	-E	output undefined voxels as 1.0e-37 (default 0.0)
            % N.B.:	default output filename = <image>_sqrt	

            assert(ischar(args));
            [s,r] = dbbash(sprintf('sqrt_4dfp %s', args));
        end
        function [s,r] = sif_4dfp__(~, args)
            %% SIF_4DFP__
            % $Id: sif_4dfp,v 1.0 Tue Feb 25 11:27:39 CST 2014 suy $
            % Usage: sif_4dfp sifstr outroot

            assert(ischar(args));
            [s,r] = dbbash(sprintf('sif_4dfp1 %s', args));
        end
        function [s,r] = t4_ident__(~, args)
            %% T4_IDENT__
            % $Id: t4_ident.c,v 1.1 2007/05/01 02:00:12 avi Exp $ 
            % Usage:  t4_ident <t4file> 
            % e.g.,   t4_ident vm11b_mpr1_to_711-2B_t4
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_ident %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.t4_ident__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = t4_inv__(~, args)
            %% T4_INV__
            % $Id: t4_inv.c,v 1.1 2007/05/01 01:20:20 avi Exp $
            % Usage:	t4_inv <t4file> [inv_t4file]
            % e.g.,	    t4_inv vm11b_anat_ave_to_vm11b_234-3_t4 [vm11b_234-3_to_vm11b_anat_ave_t4]
            % 	        option
            % 	        -u	suppress (intensity) scale field in output
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_inv %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.t4_inv__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = t4_mul__(~, args)
            %% T4_MUL__
            % $Id: t4_mul.c,v 1.1 2007/04/30 17:37:57 avi Exp $
            % Usage:	t4_mul <left_t4file> <right_t4file> [product_t4file]
            % e.g.:	    t4_mul vm11b_anat_ave_to_vm11b_234-3_t4 vm11b_234-3_to_711-2B_t4 [vm11b_anat_ave_to_711-2B_t4]
            % 	        option
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_mul %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.t4_mul__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = t4img_4dfp__(~, args)
            %% T4IMG_4DFP__
            % $Id: t4img_4dfp,v 1.14 2015/11/03 04:45:16 avi Exp $
            % Usage:	t4img_4dfp <t4file> <imgfile> [outfile]
            %  e.g.,	t4img_4dfp  vce1_mprS_to_711-2B_t4	vce1_mprS.4dfp.img -O222
            %    or,	t4img_4dfp  vce1_mprS_to_711-2B_t4 	vce1_mprS vce_mprS_711-2B -O222
            %    or,	t4img_4dfp  none			vce1_mprS vce1_mprS_222 -O222
            % N.B.:	    4dfp filename extensions are optional
            % N.B.:	    default output filename = <imgfile>t
            % N.B.:	    t4img_4dfp is a wrapper for t4imgs_4dfp; options are listed by t4imgs_4dfp usage
            % N.B.:	    option -n causes fidl ROI names to be copied to the output ifh
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4img_4dfp %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.t4img_4dfp__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = t4imgs_4dfp__(~, args)
            %% T4IMGS_4DFP__
            % $Id: t4imgs_4dfp.c,v 1.37 2015/11/03 03:50:32 avi Exp $
            % Usage:	t4imgs_4dfp [options] <inlist> <outfile>
            % 	option
            % 	-z	normalize by sqrt(n) rather than n (for z images)
            % 	-s	interpolate by 3D cubic spline (default is 3D linear)
            % 	-N	output NaN (default 0.0) for undefined values
            % 	-B	internally convert to_711-2A_t4->to_711-2B_t4
            % 	-n	use nearest neighbor interpolation
            % 	-R	suppress creation of rec file
            % 	-O111	output in 111 space instead of default 333.0 space
            % 	-O222	output in 222 space instead of default 333.0 space
            % 	-O333.n	output in 333.n space (y shifted up by n pixels)
            % 	-Omy_image	duplicate dimensions of my_image.4dfp.ifh
            % 	-@<b|l>	output big or little endian (default CPU endian)
            % N.B.:	t4file intensity scale ingnored with option -n
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4imgs_4dfp %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.t4imgs_4dfp__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = t4_resolve__(~, args)
            %% T4_RESOLVE__
            % $Id: t4_resolve.c,v 1.6 2013/09/12 00:31:12 avi Exp $
            % n=0
            % Usage:	t4_resolve <image1> <image2> ...
            % 	        option
            % 	        -v	verbose mode
            % 	        -m	generate mat file output
            % 	        -s	include intensity scale factor in t4 file output
            % 	        -w	weight inversely in proportion to scale in sub file output (sum counts mode)
            % 	        -o<str>	write resolved output with specified fileroot
            % 	        -r<flt>	set VOI rms radius in mm (default=50)
            % N.B.:	    t4_resolve looks for t4 files <image1>_to_<image2>_t4, <image1>_to_<image3>_t4, ...
            % N.B.:	    t4_resolve automatically strips filename extensions when constructing t4 filenames
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('t4_resolve %s', args));
            if (0 ~= s)
                error('mlfourdfp:abnormalExit', ...
                    'FourdfpVisitor.t4_resolve__:\n    s->%i\n    r->%s', s, r);
            end
        end
        function [s,r] = zero_slice_4dfp__(~, args)
            %% ZERO_SLICE_4DFP__
            % $Id: zero_slice_4dfp.c,v 1.2 2010/08/20 23:21:29 avi Exp $
            % Usage:	zero_slice_4dfp <(4dfp) image>
            %  e.g.,	zero_slice_4dfp vce20_mpr -z1to3
            %    or,	zero_slice_4dfp vce20_mpr <x|y|z> istart iend [outroot]
            % 	option
            % 	-<x|y|z><int>to<int>	specify x y z limits (single required argument mode)
            % 	-f	                    interpret slice numbers using 4dfp<->analyze flips
            % 	-o	                    specify output fileroot (default = <image>z)
            % 	-@<b|l>	                output big or little endian (default input endian)
            % N.B.:	slices count from 1
            % N.B.:	two usages are supported: 1 or 4 required arguments
            
            assert(ischar(args));
            [s,r] = dbbash(sprintf('zero_slice_4dfp %s', args));
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

