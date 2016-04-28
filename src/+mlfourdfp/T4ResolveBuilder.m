classdef T4ResolveBuilder 
	%% T4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Mar-2016 21:29:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	
    
    properties
        atlasTag = 'TRIO_Y_NDC'
        firstCrop = 0.5
        initialT4 = fullfile(getenv('RELEASE'), 'T_t4') % fullfile(getenv('RELEASE'), 'S_t4')
    end

	properties (Dependent)
        buildVisitor
 		product
        sessionData
        referenceImage
        referenceWeight
        sourceImage
        sourceWeight
        xfm
        
        petBlur
    end
    
    methods %% GET/SET
        function this = set.buildVisitor(this, v)
            assert(isa(v, 'mlfourdfp.FourdfpVisitor'));
            this.buildVisitor_ = v;
        end
        function v    = get.buildVisitor(this)
            v = this.buildVisitor_;
        end
        function this = set.product(this, p)
            if (isa(p, 'mlfourdfp.ImagingContext')) % preserves MR/PETImagingContexts
                this.product_ = p;
                return
            end
            this.product_ = mlfourdfp.ImagingContext(p);
        end
        function prod = get.product(this)
            prod = this.product_;
        end
        function this = set.referenceWeight(this, w)
            if (isa(w, 'mlfourdfp.ImagingContext')) % preserves MR/PETImagingContexts
                this.referenceWeight_ = w;
                return
            end
            this.referenceWeight_ = mlfourdfp.ImagingContext(w);
        end
        function w    = get.referenceWeight(this)
            % may be empty
            w = this.referenceWeight_;
        end
        function this = set.referenceImage(this, ref)
            if (isa(ref, 'mlfourdfp.ImagingContext')) % preserves MR/PETImagingContexts
                this.referenceImage_ = ref;
                return
            end
            this.referenceImage_ = mlfourdfp.ImagingContext(ref);
        end
        function ref  = get.referenceImage(this)
            ref = this.referenceImage_;
        end
        function sd   = get.sessionData(this)
            sd = this.sessionData_;
        end
        function this = set.sourceImage(this, si)
            if (isa(si, 'mlfourdfp.ImagingContext')) % preserves MR/PETImagingContexts
                this.sourceImage_ = si;
                return
            end
            this.sourceImage_ = mlfourdfp.ImagingContext(si);
        end
        function si = get.sourceImage(this)
            si = this.sourceImage_;
        end
        function this = set.sourceWeight(this, sw)
            if (isa(sw, 'mlfourdfp.ImagingContext')) % preserves MR/PETImagingContexts
                this.sourceWeight_ = sw;
                return
            end
            this.sourceWeight_ = mlfourdfp.ImagingContext(sw);
        end
        function sw = get.sourceWeight(this)
            sw = this.sourceWeight_;
        end        
        function this = set.xfm(this, x)
            %% SET.XFM casts its argument to a f. q. filename ending in this.buildVisitor.XFM_SUFFIX
            
            this.xfm_ = [myfileprefix(x) this.buildVisitor.XFM_SUFFIX];
        end
        function x    = get.xfm(this)
            %% GET.XFM 
            %  @return:
            %  - f.q. filename obtained from this.product, ending in this.buildVisitor.XFM_SUFFIX
            %  - f.q. filename with FnirtVisitor.WARPCOEF_SUFFIX replaced by XFM_SUFFIX

            import mlfsl.*;
            if (isempty(this.xfm_))
                this.xfm_ = this.product.fqfileprefix; 
            end
            lx = length(FlirtVisitor.XFM_SUFFIX);
            if (~strfind(this.xfm_(end-lx+1:end), this.buildVisitor.XFM_SUFFIX));
                this.xfm_ = [this.xfm_            this.buildVisitor.XFM_SUFFIX]; 
            end
            x = this.xfm_;
        end
        function pb   = get.petBlur(this)
            pb = this.sessionData.petBlur;
        end
    end

    methods (Abstract)
        t4ResolveSubject(this)
        t4ResolvePET(this)
    end
    
	methods
        function fn   = frameRegLoggerFilename(this, tag)
            fn = fullfile(this.sessionData_.sessionPath, [tag '_T4ResolveBuilder_frameReg_Logger_' datestr(now, 30) '.log']);
        end
        function fn   = t4ResolveAndPasteLoggerFilename(this, tag)
            fn = fullfile(this.sessionData_.sessionPath, [tag '_T4ResolveBuilder_t4ResolveAndPaste_Logger_' datestr(now, 30) '.log']);
        end        
        
        function this = t4ResolveMR(this)
        end
        function this = t4ResolveMultispectral(this)
        end
        function this = t4ResolveIterative(this, fdfp0, fdfp1, mpr)
            fdfp0 = mybasename(fdfp0);
            fdfp1 = mybasename(fdfp1);
            frame0 = 4;
            frameF = this.readFrameEnd(fdfp0);
            this.t4ResolveIteration( ...
                fdfp0, fdfp1, mpr, ...
                'frame0', frame0, 'crop', this.firstCrop);
            this.t4ResolveIteration( ...
                sprintf('%s_frames%ito%i_resolved', fdfp1, frame0, frameF), ...
                sprintf('%sr1', fdfp1), ...
                mpr, ...
                'frame0', 1, 'frameF', frameF-frame0+1, 'crop', 1);
            this.t4ResolveIteration( ...
                sprintf('%sr1_frames%ito%i_resolved', fdfp1, 1, frameF-frame0+1), ...
                sprintf('%sr2', fdfp1), ...
                mpr, ...
                'frame0', 1, 'frameF', frameF-frame0+1, 'crop', 1);
            
            this.product_ = [this.product_ sprintf('%sr2.4dfp.img', fdfp1)];
        end
        function f    = readFrameEnd(~, fdfp)
            [~,f] = mlbash(sprintf('awk ''/matrix size \\[4\\]/{print $NF}'' %s.4dfp.ifh', fdfp));
            f = str2double(f);
        end
        function f    = readSize(~, varargin)
            ip = inputParser;
            addRequired(ip, 'fdfp',  @ischar);
            addOptional(ip, 'mu', 4, @isnumeric);
            parse(ip, varargin{:});
            [~,f] = mlbash(sprintf('awk ''/matrix size \\[%i\\]/{print $NF}'' %s.4dfp.ifh', ip.Results.mu, ip.Results.fdfp));
            f = str2double(f);
        end
        function this = t4ResolveIteration(this, varargin)
            ip = inputParser;
            addRequired( ip, 'fdfp0',     @(x) lexist([x '.4dfp.img'], 'file'));
            addRequired( ip, 'fdfp1',     @ischar);
            addRequired( ip, 'mprage',    @(x) lexist([x '.4dfp.img']));
            addParameter(ip, 'frame0', 4, @isnumeric);
            addParameter(ip, 'frameF', ...
                this.readFrameEnd(varargin{1}), ...
                                          @isnumeric);
            addParameter(ip, 'crop',   1, @isnumeric);
            addParameter(ip, 'atlas', 'TRIO_Y_NDC', @(x) lexist(fullfile(getenv('REFDIR'), [x '.4dfp.img'])));
            addParameter(ip, 'blur',   5.5, @isnumeric);
            parse(ip, varargin{:});            
            
            this = this.crop(ip.Results);
            
            if (~lexist(sprintf('%s_sumt_mskt.4dfp.img', ip.Results.fdfp1), 'file'))
                this = this.msktgenInitial(ip.Results);
            end
            
            if (~lexist([ip.Results.fdfp1 '_frame' ip.Results.frameF '.4dfp.img'], 'file'))
                frameFps = this.extractFrames(ip.Results);
            end
            
            this = this.frameReg(         ip.Results, frameFps);
            this = this.t4ResolveAndPaste(ip.Results);  
            this = this.msktgenResolved(  ip.Results);
            this = this.t4ResolveTeardown(ip.Results);
        end
        function this = crop(this, t4ri)
            if (0 == t4ri.crop || t4ri.crop == 1)
                if (~lexist([t4ri.fdfp1 '.4dfp.hdr']))
                    dbbash(sprintf('copy_4dfp %s %s', t4ri.fdfp0, t4ri.fdfp1));
                end
                return
            end
                dbbash(sprintf('crop_4dfp_.sh %g %s %s', t4ri.crop, t4ri.fdfp0, t4ri.fdfp1));
        end
        function fps  = extractFrames(this, t4ri)
            frameEnd = this.readFrameEnd(t4ri.fdfp0);
            fps = cell(1, frameEnd);
            for f = t4ri.frame0:t4ri.frameF
                this.buildVisitor.extract_frame_4dfp(t4ri.fdfp1, f);
                fps{f} = sprintf('%s_frame%i', t4ri.fdfp1, f);
            end
        end
        function this = frameReg(this, t4ri, frameFps)
            
            maskfp = sprintf('%s_sumt_mskt', t4ri.fdfp1);
            file_frame0 =  sprintf('%s_frame%i.4dfp.img', t4ri.fdfp1, t4ri.frame0);
            if (~lexist(file_frame0))
                error('mlfourdfp:missingFile', 'T4ResolveBuilder.frameReg could not find %s', file_frame0);
            end     
            file_maskfp = sprintf('%s.4dfp.img', maskfp);
            if (~lexist(file_maskfp))
                error('mlfourdfp:missingFile', 'T4ResolveBuilder.frameReg could not find %s', file_maskfp);
            end
            
            if (t4ri.blur > 0)
                blurredFps = this.blurFrames(frameFps);
            end
            
            log = mlpipeline.Logger(this.frameRegLoggerFilename(t4ri.fdfp1), this);
            
            for m = t4ri.frame0:length(blurredFps)
                for n = t4ri.frame0:length(blurredFps)
                    if (m ~= n)
                        t4file = sprintf('%s_to_%s_t4', frameFps{m}, frameFps{n});
                        this.buildVisitor.imgreg_4dfp(blurredFps{m}, maskfp, blurredFps{n}, 'none', t4file, 2051, log.fqfilename);
                        this.buildVisitor.imgreg_4dfp(blurredFps{m}, maskfp, blurredFps{n}, 'none', t4file, 2051, log.fqfilename);
                    end
                end
            end
        end
        function this = t4ResolveAndPaste(this, t4ri)
            imgFns = '';
            for f = t4ri.frame0:t4ri.frameF
                imgFns = [imgFns sprintf(' %s_frame%i.4dfp.img', t4ri.fdfp1, f)];
            end
            
            log = mlpipeline.Logger(this.t4ResolveAndPasteLoggerFilename(t4ri.fdfp1), this);
            
            resolveSuff = 'resolved';
            fprintf('t4_resolve -v -m -s -o%s %s &> %s', resolveSuff, imgFns, log.fqfilename);
            this.buildVisitor.t4_resolve(resolveSuff, imgFns, 'options', '-v -m -s', 'log', log.fqfilename);
            this.t4imgFrames(t4ri, resolveSuff);
            this.pasteFrames(t4ri, resolveSuff);

            resolved = sprintf('%s_%s_%s', resolveSuff, t4ri.fdfp1, datestr(now, 30)); 
            copyfile('resolved.mat0', [resolved '.mat0']);
            copyfile('resolved.sub',  [resolved '.sub']);
        end
        function this = t4imgFrames(this, t4ri, tag)
            for f = t4ri.frame0:t4ri.frameF
                frameFile = sprintf('%s_frame%i', t4ri.fdfp1, f);
                this.buildVisitor.t4img_4dfp( ...
                    sprintf('%s_to_%s_t4', frameFile, tag), ...
                    frameFile, ...
                    sprintf('%s_%s', frameFile, tag), ...
                    'options', ['-O' frameFile]);
            end
        end
        function this = pasteFrames(this, t4ri, tag)
            
            pasteList = [t4ri.fdfp1 '_paste.lst'];
            if (lexist(pasteList)); delete(pasteList); end
            fid = fopen(pasteList, 'w');
            for f = t4ri.frame0:t4ri.frameF
                fprintf(fid, '%s_frame%i_%s.4dfp.img\n', t4ri.fdfp1, f, tag);
            end
            fclose(fid);
            
            taggedFile = sprintf('%s_frames%ito%i_%s', t4ri.fdfp1, t4ri.frame0, t4ri.frameF, tag);
            this.buildVisitor.paste_4dfp(pasteList, taggedFile, 'options', '-a ');
        end
        function fps  = blurFrames(this, fps)
            for f = 1:length(fps)
                if (~isempty(fps{f}))
                    this.buildVisitor.imgblur_4dfp(fps{f}, 5.5);
                    fps{f} = [fps{f} '_b55'];
                end
            end
        end
        function this = msktgenInitial(this, t4ri)
            %% MSKTGENINITIAL generates the 4dfp mask for dynamic imaging named
            %  [t4ri.fdfp1 '_sumt_mskt']
            %  @params t4ri is ip.Results from t4ResolveIteration            
            
            name      = t4ri.fdfp1;
            nameSumt  = [name '_sumt'];
            nameSumtG = [nameSumt '_g11'];
            mpr       = t4ri.mprage;
            atl       = t4ri.atlas;
            log       = sprintf('msktgenInitial_%s.log', datestr(now, 30));
            
            mprToAtlT4 = [mpr '_to_' atl '_t4'];
            if (~lexist(mprToAtlT4, 'file'))
                error('mlfourdfp:missingPrerequisiteFile', 'T4ResolveBuilder.msktgenInitial:  file %s not found', mprToAtlT4);
            end
            
            nameToMprT4 = [name '_to_' mpr '_t4'];
            if (lexist(nameToMprT4, 'file'))
                delete(nameToMprT4);
            end
            this.buildVisitor.t4_inv(this.initialT4, nameToMprT4);
            this.buildVisitor.actmapf_4dfp( ...
                sprintf('%i+', this.readFrameEnd(t4ri.fdfp0)), name, 'options', '-asumt');
            this.buildVisitor.gauss_4dfp(nameSumt, 1.1, nameSumtG); 
            
            petMsk = this.maskBoundaries(nameSumtG);
            
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, petMsk, nameToMprT4, 4099, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, petMsk, nameToMprT4, 4099, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, petMsk, nameToMprT4, 3075, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, petMsk, nameToMprT4, 2051, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, petMsk, nameToMprT4, 2051, log); % 10243
            this.buildVisitor.t4img_4dfp( nameToMprT4, nameSumtG, [name '_on_' mpr], 'options', ['-n -O' mpr]);            
            this.buildVisitor.t4_mul(     nameToMprT4, mprToAtlT4, [nameSumt '_to_' atl '_t4']);
            
            this.buildVisitor.msktgen_4dfp(nameSumt, 20, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
        end
        function msk  = maskBoundaries(this, fdfp)
            this.buildVisitor.nifti_4dfp_n(fdfp);
            ic  = mlfourd.ImagingContext(fdfp);
            ic  = ic.ones;
            ic.noclobber = 0;
            ic.saveas([fdfp '_ones']);
            this.buildVisitor.nifti_4dfp_4([fdfp '_ones']);
            msk = this.zeroSlicesOnBoundaries([fdfp '_ones'], 3);
        end
        function fdfp = zeroSlicesOnBoundaries(this, fdfp, n)            
            N     = this.readSize(fdfp, 3);
            fdfp1 = sprintf('%s_z%ito%i', fdfp,  1, n);
            fdfp2 = sprintf('%s_z%ito%i', fdfp1, N-n+1, N);
            this.buildVisitor.zero_slice_4dfp(fdfp,  'z', 1, n,     fdfp1);
            this.buildVisitor.zero_slice_4dfp(fdfp1, 'z', N-n+1, N, fdfp2);
            fdfp  = fdfp2;            
        end
        function this = msktgenMprage(this, fdfp, atl)
            log = sprintf('msktgenMprage_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fdfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
            this.buildVisitor.msktgen_4dfp(fdfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
        end
        function this = msktgenResolved(this, t4ri)

            this = this.msktgenInitial(t4ri);
            
            % for viewing, inspection   
            name = t4ri.fdfp1;
            mask = [name '_sumt_mskt'];
            this.buildVisitor.imgblur_4dfp( name, 5.5);
            this.buildVisitor.maskimg_4dfp([name '_b55'], mask, [name '_b55_mskt']);
        end
        function this = t4ResolveTeardown(this, t4ri)
            return
            
            name = t4ri.fdfp1;
            for f = 1:t4ri.frameF
                delete(sprintf('%s_frame%i.4dfp.*', name, f));
            end
            for f = t4ri.frame0:t4ri.frameF
                delete(sprintf('%s_frame%i_b55.4dfp.*', name, f));
                delete(sprintf('%s_frame%i_resolved.4dfp.*', name, f));
            end
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
        
 		function this = T4ResolveBuilder(varargin)
 			%% T4RESOLVEBUILDER
 			%  Usage:  this = T4ResolveBuilder()
 			
            if (0 == nargin); return; end
            
            %% invoke copy-ctor
            
            if (1 == nargin && isa(varargin{1}, 'mlfourdfp.T4ResolveBuilder'))
                arg = varargin{1};                
                this.sessionData_  = arg.sessionData_;                
                this.buildVisitor_ = arg.buildVisitor_;
                return
            end
            
            %% manage parameters 
            
            import mlfourdfp.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'sessionData',  [],             @(x) isa(x, 'mlpipeline.ISessionData'));
            addParameter(ip, 'buildVisitor', FourdfpVisitor, @(x) isa(x, 'mlfourdfp.FourdfpVisitor'));
            parse(ip, varargin{:});
            
            this.sessionData_  = ip.Results.sessionData;
            this.buildVisitor_ = ip.Results.buildVisitor;
 		end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)        
        buildVisitor_
        product_
        referenceImage_
        referenceWeight_
        sessionData_
        sourceImage_
        sourceWeight_
        xfm_
    end
    
    methods (Access = protected)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

 