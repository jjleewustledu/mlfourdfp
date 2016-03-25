classdef T4ResolveBuilder 
	%% T4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Mar-2016 21:29:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	
    
    properties
        blur = 0
    end

	properties (Dependent)
        buildVisitor
        frameRegLoggerFilename
        t4ResolveAndPasteLoggerFilename
 		product
        sessionData
        referenceImage
        referenceWeight
        sourceImage
        sourceWeight
        xfm
    end
    
    methods %% GET/SET
        function this = set.buildVisitor(this, v)
            assert(isa(v, 'mlfourdfp.FourdfpVisitor'));
            this.buildVisitor_ = v;
        end
        function v    = get.buildVisitor(this)
            v = this.buildVisitor_;
        end
        function fn   = get.frameRegLoggerFilename(this)
            fn = fullfile(this.sessionData_.sessionPath, ['T4ResolveBuilder_frameReg_Logger_' datestr(now, 30) '.log']);
        end
        function fn   = get.t4ResolveAndPasteLoggerFilename(this)
            fn = fullfile(this.sessionData_.sessionPath, ['T4ResolveBuilder_t4ResolveAndPaste_Logger_' datestr(now, 30) '.log']);
        end
        function this = set.product(this, p)
            if (isa(p, 'mlfourd.ImagingContext')) % preserves MR/PETImagingContexts
                this.product_ = p;
                return
            end
            this.product_ = mlfourd.ImagingContext(p);
        end
        function prod = get.product(this)
            prod = this.product_;
        end
        function this = set.referenceWeight(this, w)
            if (isa(w, 'mlfourd.ImagingContext')) % preserves MR/PETImagingContexts
                this.referenceWeight_ = w;
                return
            end
            this.referenceWeight_ = mlfourd.ImagingContext(w);
        end
        function w    = get.referenceWeight(this)
            % may be empty
            w = this.referenceWeight_;
        end
        function this = set.referenceImage(this, ref)
            if (isa(ref, 'mlfourd.ImagingContext')) % preserves MR/PETImagingContexts
                this.referenceImage_ = ref;
                return
            end
            this.referenceImage_ = mlfourd.ImagingContext(ref);
        end
        function ref  = get.referenceImage(this)
            ref = this.referenceImage_;
        end
        function sd   = get.sessionData(this)
            sd = this.sessionData_;
        end
        function this = set.sourceImage(this, si)
            if (isa(si, 'mlfourd.ImagingContext')) % preserves MR/PETImagingContexts
                this.sourceImage_ = si;
                return
            end
            this.sourceImage_ = mlfourd.ImagingContext(si);
        end
        function si = get.sourceImage(this)
            si = this.sourceImage_;
        end
        function this = set.sourceWeight(this, sw)
            if (isa(sw, 'mlfourd.ImagingContext')) % preserves MR/PETImagingContexts
                this.sourceWeight_ = sw;
                return
            end
            this.sourceWeight_ = mlfourd.ImagingContext(sw);
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
    end

	methods
        function this = t4ResolveSubject(this)
            this.sourceImage = this.sessionData.pet;
            this = this.t4ResolvePET;
            pet  = this.product;

            this.sourceImage = this.sessionData.mr;
            this = this.t4ResolveMR;
            mr   = this.product;            

            this.sourceImage = mlfourd.ImagingContext({pet mr});
            this = this.t4Resolve;
        end
        function this = t4ResolvePET(this)

            atlTag  = 'TRIO_Y_NDC';
            subject = 'NP995_09';
            convdir = fullfile(getenv('PPG'), 'converted', subject, '');
            nacdir  = fullfile(convdir, 'NAC', '');
            workdir = fullfile(getenv(PPG), 'jjlee', subject, 'NAC', '');
            mpr     = [subject '_mpr'];

            if (~lexist([mpr '.4dfp.img'], 'file'))
                this.buildVisitor.lns_4dfp(fullfile(convdir, mpr));
            end
            if (~lexist([mpr '_to_' atlTag '_t4']))
                this.msktgenMprage(mpr, atlTag);
            end

            tracer = 'fdg';
            for visit = 1:2
                tracerdir = fullfile(workdir, sprintf('%s_v%i', upper(tracer), visit));
                this.buildVisitor.mkdir(tracerdir);
                this.buildVisitor.pushd(tracerdir);
                this.buildVisitor.lns(     fullfile(workdir, [mpr '_to_' atlTag '_t4']));
                this.buildVisitor.lns_4dfp(fullfile(convdir,  mpr));
                this.buildVisitor.lns_4dfp(fullfile(nacdir, sprintf('%s%s_v%i_NAC', subject, upper(tracer), visit)));
                fdfp0 = sprintf('%s%s_v%i_NAC', subject, upper(tracer), visit);
                fdfp1 = sprintf('%sv%i', tracer, visit);
                dbbash( sprintf('t4_resolver_top_.sh %s %s &', fdfp0, fdfp1));
                this.buildVisitor.popd;
            end
            
            tracers = {'ho' 'oo'};
            for t = 1:length(tracers)
                for visit = 1:2
                    for scan = 1:2
                        tracerdir = fullfile(workdir, sprintf('%s%i_v%i', upper(tracer), scan, visit));
                        this.buildVisitor.mkdir(tracerdir);
                        this.buildVisitor.pushd(tracerdir);
                        this.buildVisitor.lns(     fullfile(workdir, [mpr '_to_' atlTag '_t4']));
                        this.buildVisitor.lns_4dfp(fullfile(convdir,  mpr));
                        this.buildVisitor.lns_4dfp(fullfile(nacdir, sprintf('%s%s%i_v%i_NAC', subject, upper(tracer), scan, visit)));
                        fdfp0 = sprintf('%s%s%i_v%i_NAC', subject, upper(tracer), scan, visit);
                        fdfp1 = sprintf('%s%iv%i', tracer, scan, visit);
                        dbbash( sprintf('t4_resolver_top_.sh %s %s &', fdfp0, fdfp1));
                        this.buildVisitor.popd;
                    end
                end
            end
        end
        function this = t4ResolveMR(this)
        end
        function this = t4ResolveIterative(this, fdfp0, fdfp1)
            fdfp0 = basename(fdfp0);
            fdfp1 = basename(fdfp1);
            frame0 = 4;
            frameF = this.readFrameEnd(fdfp0);
            this.t4ResolveIteration( ...
                fdfp0, fdfp1, ...
                'frame0', frame0, 'crop', 0.5);
            this.t4ResolveIteration( ...
                sprintf('%s_frames%sto%s_resolved', fdfp1, frame0, frameF), ...
                sprintf('%sr1', fdfp1), ...
                'frame0', frame0, 'frameF', frameF,'crop', 1);
            this.t4ResolveIteration( ...
                sprintf('%sr1_frames%sto%s_resolved', fdfp1, 1, frameF-frame0+1), ...
                sprintf('%sr2', fdfp1), ...
                'frame0', 1, 'frameF', frameF-frame0+1, 'crop', 1);
        end
        function f    = readFrameEnd(~, fdfp)
            [~,f] = mlbash(sprintf('awk ''/matrix size \\[4\\]/{print $NF}'' %s.4dfp.ifh', fdfp));
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
            addParameter(ip, 'crop',   1, @islogical);
            addParameter(ip, 'atlas', 'TRIO_Y_NDC', @(x) lexist(fullfile(getenv('REFDIR'), [x '.4dfp.img'])));
            addParameter(ip, 'blur',   5.5, @isnumeric);
            parse(ip, varargin{:});            
            
            this = this.crop(ip.Results);
            this = this.msktgenInitial(ip.Results);
            
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
                dbbash(sprintf('copy_4dfp %s %s', t4ri.fdfp0, t4ri.fdfp1));  
                return
            end
                dbbash(sprintf('crop_4dfp_.sh %s %s %s', t4ri.crop, t4ri.fdfp0, t4ri.fdfp1));
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
            file_frame0 =  sprintf('%s_frame%i.4dfp.img', t4ri.fdfp0, t4ri.frame0);
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
            
            log = mlpipeline.Logger(this.frameRegLoggerFilename, this);
            
            for m = 1:length(blurredFps)
                for n = 1:length(blurredFps)
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
                imgFns = [imgFns sprintf('%s_frame%i.4dfp.img', t4ri.fdfp1, f)];
            end
            
            log = mlpipeline.Logger(this.t4ResolveAndPasteLoggerFilename, this);
            
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
            frameN = this.readFrameEnd(t4ri.fdfp0);
            for f = t4ri.frame0:frameN
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
                fprintf(fid, '%s_frame%i_%s.4dfp.img', t4ri.fdfp1, f, tag);
            end
            fclose(fid);
            
            taggedFile = sprintf('%s_frames%ito%i_%s', t4ri.fdfp1, t4ri.frame0, t4ri.frameF, tag);
            this.buildVisitor.paste_4dfp(pasteList, taggedFile, 'options', '-a ');
        end
        function fps  = blurFrames(this, fps)
            for f = 1:length(fps)
                this.buildVisitor.imgblur_4dfp(fps{f}, 5.5);
                fps{f} = [fps{f} '_b55'];
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
            
            if (lexist(sprintf('%s_mskt.4dfp.img', nameSumt), 'file'))
                warning('mlfourdfp:nothingToDo', 'T4ResolveBuilder.msktgenInitial');
                return
            end
            mprToAtlT4 = [mpr '_to_' atl '_t4'];
            if (~lexist(mprToAtlT4, 'file'))
                error('mlfourdfp:missingPrerequisiteFile', 'T4ResolveBuilder.msktgenInitial:  file %s not found', mprToAtlT4);
            end
            
            nameToMprT4 = [name '_to_' mpr '_t4'];
            this.buildVisitor.t4_inv(fullfile(getenv('RELEASE'), 'S_t4'), nameToMprT4);
            this.buildVisitor.actmapf_4dfp( ...
                sprintf('%i+', this.readFrameEnd(t4ri.fdfp0)), name, 'options', '-asumt');
            this.buildVisitor.gauss_4dfp(nameSumt, 1.1, nameSumtG); 
            
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, 'none', nameToMprT4, 4099, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, 'none', nameToMprT4, 4099, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, 'none', nameToMprT4, 3075, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, 'none', nameToMprT4, 2051, log);
            this.buildVisitor.imgreg_4dfp(mpr, 'none', nameSumtG, 'none', nameToMprT4, 2051, log);
            this.buildVisitor.t4img_4dfp( nameToMprT4, nameSumtG, [name '_on_' mpr], 'options', ['-n -O' mpr]);            
            this.buildVisitor.t4_mul(     nameToMprT4, mprToAtlT4, [name '_to_' atl '_t4']);
            
            this.buildVisitor.imgblur_4dfp( nameSumt, 5.5);
            this.buildVisitor.msktgen_4dfp([nameSumt '_b55'], 20, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
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
            this.buildVisitor.maskimg_4dfp([name '_b55'], mask, [name '_b55_sumt']);
        end
        function this = t4ResolveTeardown(this, t4ri)
            
            name = t4ri.fdfp1;
            for f = 1:t4ri.frameF
                delete(sprintf('%s_frame%i.4dfp.*', name, f));
            end
            for f = t4ri.frame0:t4ri.frameF
                delete(sprintf('%s_frame%i_b55.4dfp.*', name, f));
                delete(sprintf('%s_frame%i_resolved.4dfp.*', name, f));
            end
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

