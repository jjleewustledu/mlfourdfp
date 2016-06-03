classdef T4ResolveBuilder 
	%% T4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Mar-2016 21:29:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	
    
    properties (Abstract)
        atlasTag
        blurArg
        firstCrop
        gaussArg
        initialT4 
    end

    methods (Abstract)
        t4ResolveSubject(this)
        t4ResolvePET(this)
    end
    
    properties (Constant)
        F_HALF_x_FWHM = 10*0.4412712
    end
    
    properties
        framesToSkip = 3
        keepForensics = false
    end
    
	properties (Dependent)
        buildVisitor
        petBlur
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
        function si   = get.sourceImage(this)
            si = this.sourceImage_;
        end
        function this = set.sourceWeight(this, sw)
            if (isa(sw, 'mlfourdfp.ImagingContext')) % preserves MR/PETImagingContexts
                this.sourceWeight_ = sw;
                return
            end
            this.sourceWeight_ = mlfourdfp.ImagingContext(sw);
        end
        function sw   = get.sourceWeight(this)
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
    
	methods
        function this = t4ResolveMR(this)
        end
        function this = t4ResolveMultispectral(this)
        end
        function this = t4ResolveIterative(this, fdfp0, fdfp1, mpr)
            fdfp0  = mybasename(fdfp0);
            fdfp1  = mybasename(fdfp1);
            frame0 = this.framesToSkip+1;
            frameF = this.readFrameEnd(fdfp0);
            this = this.t4ResolveIterate( ...
                fdfp0, fdfp1, ...
                'mprage', mpr, ...
                'frame0', frame0, 'crop', this.firstCrop);
            this = this.t4ResolveIterate( ...
                sprintf('%s_frames%ito%i_resolved', fdfp1, frame0, frameF), ...
                sprintf('%sr1', fdfp1), ...
                'mprage', mpr, ...
                'frame0', 1, 'frameF', frameF-frame0+1);
            
            return
            
            this = this.t4ResolveIterate( ...
                sprintf('%sr1_frames%ito%i_resolved', fdfp1, 1, frameF-frame0+1), ...
                sprintf('%sr2', fdfp1), ...
                'mprage', mpr, ...
                'frame0', 1, 'frameF', frameF-frame0+1);             
            this = this.t4ResolveFinalize( ...
                sprintf('%sr2_frames%ito%i_resolved', fdfp1, 1, frameF-frame0+1), ...
                sprintf('%s_resolvedFinal', fdfp1), ...
                'mprage', mpr, ...
                'frame0', 1, 'frameF', frameF-frame0+1);
            this.product_ = [this.product_ this.filename([fdfp1 'r2'])];
        end
        function this = t4ResolveIterate(this, varargin)
            ip = this.t4ResolveParser(varargin{:});          
            fdfp1 = ip.Results.fdfp1;
            
            if (~lexist(this.filename(fdfp1)))
                this = this.crop(ip.Results);
            end
            if (~lexist(this.filename(this.fdfpMskt(this.fdfpGauss(this.fdfpSumt(fdfp1)))), 'file'))
                this = this.msktgenInitial(ip.Results);
            end
            frameFps = this.extractFrames(ip.Results);
            this = this.frameReg(ip.Results, frameFps);
            this = this.t4ResolveAndPaste(ip.Results); 
            this = this.teardownT4ResolveIteration(ip.Results);
            this.buildVisitor.imgblur_4dfp(fdfp1, this.blurArg);
        end
        function this = t4ResolveFinalize(this, varargin)
            ip = this.t4ResolveParser(varargin{:});             
            fdfp1 = ip.Results.fdfp1;       
            
            this = this.crop(ip.Results);
            if (~lexist( ...
                    this.filename(this.fdfpMskt(this.fdfpGauss(this.fdfpSumt(fdfp1)))), 'file'))
                this = this.msktgenInitial(ip.Results);
            end            
            this.buildVisitor.imgblur_4dfp(fdfp1, this.blurArg);
        end
        function fdfps  = extractFrames(this, ipr)
            frameEnd = this.readFrameEnd(ipr.fdfp0);
            fdfps = cell(1, frameEnd);
            for f = ipr.frame0:ipr.frameF
                fdfps{f} = sprintf('%s_frame%i', ipr.fdfp1, f);
                if (~lexist(this.filename(fdfps{f})))
                    this.buildVisitor.extract_frame_4dfp(ipr.fdfp1, f);
                end
            end
        end
        function this = frameReg(this, ipr, frameFdfps)
            frame0Fdfp = this.filename(sprintf('%s_frame%i', ipr.fdfp1, ipr.frame0));
            if (~lexist(frame0Fdfp))
                error('mlfourdfp:missingFile', 'T4ResolveBuilder.frameReg could not find %s', frame0Fdfp);
            end            
            blurredFdfps = this.blurFrames(ipr, frameFdfps);
            log = mlpipeline.Logger(this.loggerFilename('frameReg', ipr.fdfp1), this);            
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
        function fdfp = frameRegMask(this, ipr)
            fdfp = 'none';
            
            %fdfp = this.fdfpMskt(this.fdfpGauss(this.fdfpSumt(ipr.fdfp1)));
            %if (~lexist(this.filename(fdfp)))
            %    error('mlfourdfp:missingFile', 'T4ResolveBuilder.frameReg could not find %s', maskfn);
            %end
        end
        function this = t4ResolveAndPaste(this, ipr)
            imgFns = '';
            for f = ipr.frame0:ipr.frameF
                imgFns = [imgFns this.filename(sprintf(' %s_frame%i', ipr.fdfp1, f))];
            end
            
            log = mlpipeline.Logger(this.loggerFilename('t4ResolveAndPaste', ipr.fdfp1), this);
            
            resolveSuff = 'resolved';
            fprintf('t4_resolve -v -m -s -o%s %s &> %s', resolveSuff, imgFns, log.fqfilename);
            this.buildVisitor.t4_resolve(resolveSuff, imgFns, 'options', '-v -m -s', 'log', log.fqfilename);
            this.t4imgFrames(ipr, resolveSuff);
            this.pasteFrames(ipr, resolveSuff);

            resolved = sprintf('%s_%s_%s', resolveSuff, ipr.fdfp1, datestr(now, 30)); 
            copyfile('resolved.mat0', [resolved '.mat0']);
            copyfile('resolved.sub',  [resolved '.sub']);
        end
        function this = t4imgFrames(this, ipr, tag)
            for f = ipr.frame0:ipr.frameF
                frameFile = sprintf('%s_frame%i', ipr.fdfp1, f);
                this.buildVisitor.t4img_4dfp( ...
                    sprintf('%s_to_%s_t4', frameFile, tag), ...
                    frameFile, ...
                    'out', sprintf('%s_%s', frameFile, tag), ...
                    'options', ['-O' frameFile]);
            end
        end
        function this = pasteFrames(this, ipr, tag)
            
            pasteList = [ipr.fdfp1 '_paste.lst'];
            if (lexist(pasteList)); delete(pasteList); end
            fid = fopen(pasteList, 'w');
            for f = ipr.frame0:ipr.frameF
                fprintf(fid, '%s_frame%i_%s.4dfp.img\n', ipr.fdfp1, f, tag);
            end
            fclose(fid);
            
            taggedFile = sprintf('%s_frames%ito%i_%s', ipr.fdfp1, ipr.frame0, ipr.frameF, tag);
            this.buildVisitor.paste_4dfp(pasteList, taggedFile, 'options', '-a ');
        end
        function this = msktgenMprage(this, fdfp, atl)
            log = sprintf('msktgenMprage_%s.log', datestr(now, 30));
            this.buildVisitor.mpr2atl_4dfp(fdfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
            this.buildVisitor.msktgen_4dfp(fdfp, 'options', ['-T' fullfile(getenv('REFDIR'), atl)], 'log', log);
        end
        function this = msktgenInitial(this, ipr)
            %% MSKTGENINITIAL generates:
            %  - blurred MPRAGE named [ipr.mprage '_gNN']
            %  -         summed dynamic imaging named [ipr.fdfp1 '_sumt']
            %  - blurred summed dynamic imaging named [ipr.fdfp1 '_sumt_gNN']
            %  - t4 files for [ipr.fdfp1 '_sumt_gNN'] to ipr.mprage; 
            %                  ipr.mprage to ipr.atlas; 
            %                 [ipr.fdfp1 '_sumt_gNN'] to ipr.atlas;
            %                  ipr.fdfp1 to ipr.atlas
            %  - blurred summed dynamic imaging named [ipr.fdfp1 '_sumt_gNN_on_' ipr.mprage]
            %  -          dynamic imaging named [ipr.fdfp1 '_on_' ipr.atlas]
            %  - mask for dynamic imaging named [ipr.fdfp1 '_sumt_gNN_mskt']
            %  @params ipr is ip.Results from this.t4ResolveParser               
            
            [sumtG,toMprGT4,log] = this.alignDynamicToMprage(ipr);             
            this                 = this.alignDynamicToAtlas(ipr, sumtG, toMprGT4);    
            
            this.buildVisitor.msktgen_4dfp(sumtG, 20, 'options', ['-T' fullfile(getenv('REFDIR'), ipr.atlas)], 'log', log);
        end
        function [sumtG,toMprGT4,log] = ...
                        alignDynamicToMprage(this, ipr)
            
            mprG = this.gaussMprage(ipr);
            
            [sumt,sumtG,toMprGT4] = this.sumtDynamic(ipr, mprG);
            
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
        function this = alignDynamicToAtlas(this, ipr, sumtG, toMprGT4)
            mprToAtlT4 = [ipr.mprage '_to_' ipr.atlas '_t4'];
            if (~lexist(mprToAtlT4, 'file'))
                this = this.msktgenMprage(ipr.mprage, ipr.atlas);
            end
            dynToAtlT4 = [sumtG '_to_' ipr.atlas '_t4'];
            this.buildVisitor.t4_mul(toMprGT4, mprToAtlT4, dynToAtlT4);
            this.buildVisitor.t4img_4dfp(dynToAtlT4, ipr.fdfp1, ...
                'out', [ipr.fdfp1, '_on_' ipr.atlas], ...
                'options', ['-O' fullfile(getenv('REFDIR'), ipr.atlas)]);       
        end
        function this = teardownT4ResolveIteration(this, ipr)
            if (this.keepForensics); return; end
            
            name = ipr.fdfp1;
            for f = 1:ipr.frameF
                delete(sprintf('%s_frame%i.4dfp.*', name, f));
            end
            for f = ipr.frame0:ipr.frameF
                delete(sprintf('%s_frame%i_b55.4dfp.*', name, f));
                delete(sprintf('%s_frame%i_resolved.4dfp.*', name, f));
            end
        end
        function fqfp = buildUmapOnSumt(this)
            ct   = [this.sessionData.ct_fqfp '_on_' mybasename(this.sessionData.mpr_fqfp)];
            mpr  = this.sessionData.mpr_fqfp;
            sumt = this.fdfpSumt(this.sessionData.fdg_fqfp);
            t40  = this.buildVisitor.t4_inv(fullfile(getenv('RELEASE'), 'S_t4'), 'out', 'S_inv_t4');            
            %t4cm = this.buildVisitor.align_multiSpectral(ct,   mpr, 't40', t40);
            t4sm = this.buildVisitor.align_multiSpectral(sumt, mpr, 't40', t40);
            t4ms = this.buildVisitor.t4_inv(t4sm);
            %t4   = this.buildVisitor.t4_mul(t4cm, t4ms);
            ct   = this.buildVisitor.t4img_4dfp(t4ms, ct, 'options', ['-O' sumt]);
            fqfp = this.buildUmapFromCt(ct, this.sessionData.umap_fqfp);
        end
        function this = buildCt(this)
            % ct = ImagingContext('ct.nii')
            % ct = ct.thresh(900)
            % ct.size
            % ct.csize
            % ct.niftid.size
            % ct = ct.blurred([51 51 7])
            % ct.view
            % ct = ct.thresh(50)
            % ct = ct.binarized
            % ct.view
            % ct = ImagingContext('ct.nii')
            % ct = ct.thresh(900)
            % ct = ct.blurred([51 51 51])
            % ct.view
            % ct = ct.thresh(200)
            % ct = ct.binarized
            % ct.view('ct.nii')
            % ct.fileprefix
            % ct.save
            % ct0 = ImagingContext('ct.nii');
            % ct0 = ct0.masked(ct.niftid);
            % ct0.fileprefix = 'ct_masked';
            % ct0.save
        end
        function umap = buildUmapFromCt(this, ct, umap)
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
        function umapsF = buildUmapFrames(this, ipr)
            
            frameIndices = ipr.frame0:ipr.frameF;
            umap0  = sprintf('%s_umap', ipr.fdfp0);
            umapsF = cellfun( ...
                @(x) sprintf('%s_umap_frame%i', ipr.fdfp1, x), ...
                num2cell(frameIndices), ...
                'UniformOutput', false);  
            this.buildVisitor.copy_4dfp(myfileprefix(this.sessionData.umap_fqfn), umap0);  
            
            t4_frames_to_r0     = this.t4FramesCells( ipr.fdfp0,       frameIndices);
            t4_r0_to_r1         = this.t4FramesCells([ipr.fdfp0 'r1'], frameIndices-ipr.frame0+1);
            t4_r1_to_r2         = this.t4FramesCells([ipr.fdfp0 'r2'], frameIndices-ipr.frame0+1);            
            t4_frames_from_r0   = this.t4sInv(t4_frames_to_r0);  
            t4_r0_from_r1       = this.t4sInv(t4_r0_to_r1); 
            t4_r1_from_r2       = this.t4sInv(t4_r1_to_r2);                                
            t4_r2_from_sumt     = this.t4Cells( ...
                                  this.imgregPET(this.sumtImage(ipr.fdfp0), this.resolvedSumtImage(ipr)), ...
                                  ipr.frameF-ipr.frame0+1);
            t4_frames_from_sumt = this.t4sMuls( ...
                t4_frames_from_r0, t4_r0_from_r1, t4_r1_from_r2, t4_r2_from_sumt);
            this.t4sImg(t4_frames_from_sumt, umap0, umapsF, umap0);            
            this.pasteFramesUmap(ipr);
        end
        function this = pasteFramesAC(this, ipr)
            pasteList = [ipr.fdfp1 '_paste.lst'];
            if (lexist(pasteList)); delete(pasteList); end
            fid = fopen(pasteList, 'w');
            for f = ipr.frame0:ipr.frameF
                fprintf(fid, '%s_frame%i_OP.4dfp.img\n', ipr.fdfp1, f);
            end
            fclose(fid);
            
            taggedFile = sprintf('%s_frames%ito%i', ipr.fdfp1, 1, ipr.frameF);
            this.buildVisitor.paste_4dfp(pasteList, taggedFile, 'options', '-a ');
        end
        function this = pasteFramesUmap(this, ipr)
            ipr.fdfp1 = [ipr.fdfp1 '_umap'];
            pasteList = [ipr.fdfp1 '_paste.lst'];
            if (lexist(pasteList)); delete(pasteList); end
            fid = fopen(pasteList, 'w');
            for f = 1:ipr.frame0-1
                fprintf(fid, '%s.4dfp.img\n', ipr.fdfp1);
            end
            for f = ipr.frame0:ipr.frameF
                fprintf(fid, '%s_frame%i.4dfp.img\n', ipr.fdfp1, f);
            end
            fclose(fid);
            
            taggedFile = sprintf('%s_frames%ito%i', ipr.fdfp1, 1, ipr.frameF);
            this.buildVisitor.paste_4dfp(pasteList, taggedFile, 'options', '-a ');
        end
        function t4   = imgregPET(this, varargin)
            %% ALIGNPET blurs, calls imgreg_4dfp with mode 2051 and writes a log.
            %  @return t4 is the t4-file for the transformation.            
            
            ip           = this.t4ResolveParser(varargin{:});            
            blurredFdfp0 = this.ensureBlurred(ip.Results.fdfp0);
            blurredFdfp1 = this.ensureBlurred(ip.Results.fdfp1);
            t4           = sprintf('%s_to_%s_t4', blurredFdfp0, blurredFdfp1);
            log          = mlpipeline.Logger(this.loggerFilename('imgregPET', ip.Results.fdfp0), this);
            
            this.buildVisitor.imgreg_4dfp(blurredFdfp0, 'none', blurredFdfp1, 'none', t4, 2051, log.fqfilename);
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
        function ip    = t4ResolveParser(this, varargin)
            ip = inputParser;
            addRequired( ip, 'fdfp0',   @(x) lexist(this.filename(x), 'file'));
            addRequired( ip, 'fdfp1',   @ischar);
            addParameter(ip, 'mprage',  this.sessionData.mpr_fqfn, ...
                                        @(x) lexist(this.filename(x), 'file'));
            addParameter(ip, 'frame0',  this.framesToSkip+1, ...
                                        @isnumeric);
            addParameter(ip, 'frameF',  this.readFrameEnd(varargin{1}), ...
                                        @isnumeric);
            addParameter(ip, 'crop',    1, ...
                                        @isnumeric);
            addParameter(ip, 'atlas',   'TRIO_Y_NDC', ...
                                        @(x) lexist(fullfile(getenv('REFDIR'), this.filename(x))));
            addParameter(ip, 'blur',    this.blurArg, ...
                                        @isnumeric);
            addParameter(ip, 'rnumber', ...
                                        0, @isnumeric );
            parse(ip, varargin{:});  
        end     
        function fqfp  = transverseMpr(this)
            fqfp = [this.sessionData.mpr_fqfp '_trans'];
            if (~lexist(this.filename(fqfp), 'file'))
                fqfp = this.buildVisitor.t4img_4dfp( ...
                    fullfile(getenv('RELEASE'), 'S_t4'), this.sessionData.mpr_fqfp, ...
                    'out', fqfp, 'options', ['-O' this.sessionData.mpr_fqfp]);
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
        function fdfps = blurFrames(this, ipr, fdfps)
            if (~(ipr.blur > 0)); return; end
            for f = 1:length(fdfps)
                if (~isempty(fdfps{f}))
                    this.buildVisitor.imgblur_4dfp(fdfps{f}, this.blurArg);
                    fdfps{f} = this.fdfpBlur(fdfps{f});
                end
            end
        end
        function c    = char2cell(~, c, sz)
            if (~iscell(c))
                c = cellfun(@(x) c, cell(sz), 'UniformOutput', false);
            end
        end
        function this = crop(this, ipr)
            if (0 == ipr.crop || ipr.crop == 1)
                if (~lexist(this.filename(ipr.fdfp1)))
                    dbbash(sprintf('copy_4dfp %s %s', ipr.fdfp0, ipr.fdfp1));
                end
                return
            end
                dbbash(sprintf('cropfrac_4dfp %g %s %s', ipr.crop, ipr.fdfp0, ipr.fdfp1));
        end
        function fdfp = ensureBlurred(this, fdfp)
            blurTag = sprintf('_b%i', floor(this.blurArg*10));
            if (~lstrfind(fdfp, blurTag))
                this.buildVisitor.imgblur_4dfp(fdfp, this.blurArg);
                fdfp = this.fdfpBlur(fdfp);
            end
        end
        function fn   = filename(~, fp)
            fn = [fp '.4dfp.img'];
        end
        function fdfp = fdfpBlur(this, fdfp)
            fdfp = sprintf('%s_b%i', fdfp, floor(this.blurArg*10));
        end
        function fdfp = fdfpGauss(this, fdfp)
            fdfp = sprintf('%s_g%i', fdfp, floor(this.gaussArg*10));
        end
        function fdfp = fdfpMskt(~, fdfp)
            fdfp = [fdfp '_mskt'];
        end
        function fdfp = fdfpSumt(~, fdfp)
            fdfp = [fdfp '_sumt'];
        end
        function mprG = gaussMprage(this, ipr)
            mprG = this.fdfpGauss(ipr.mprage);            
            this.buildVisitor.gauss_4dfp(ipr.mprage, this.gaussArg, mprG);   
        end 
        function fn   = loggerFilename(this, func, tag)
            fn = fullfile(this.sessionData_.sessionPath, ...
                 sprintf('%s_T4ResolveBuilder_%s_Logger_%s.log', tag, func, datestr(now, 30)));
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
        function fdfp = resolvedSumtImage(this, ipr)
            fdfp = this.sumtImage( ...
                sprintf('%sr%i_frames%ito%i_resolved', ipr.fdfp1, ipr.rnumber, 1, ipr.frameF-ipr.frame0+1));
        end
        function [sumt,sumtG,toMprGT4] = ...
                        sumtDynamic(this, ipr, mprG)
            sumt  = this.fdfpSumt(ipr.fdfp1);
            sumtG = this.fdfpGauss(sumt);                        
            toMprGT4 = [sumtG '_to_' mprG '_t4'];
            if (lexist(toMprGT4, 'file'))
                delete(toMprGT4);
            end
            this.buildVisitor.t4_inv(this.initialT4, 'out', toMprGT4);
            this.buildVisitor.actmapf_4dfp( ...
                sprintf('%i+', this.readFrameEnd(ipr.fdfp0)), ipr.fdfp1, 'options', '-asumt');
        end
        function fdfp = sumtImage(this, fdfp0)
            fdfp = sprintf('%s_sumt', fdfp0);
            if (~lexist(this.filename(fdfp), 'file'))
                this.buildVisitor.actmapf_4dfp( ...
                    sprintf('%i+', this.readFrameEnd(fdfp0)), fdfp0, 'options', '-asumt');
            end
        end
        function t4c  = t4Cells(~, t4, N)
            t4c = cellfun( ...
                @(x) t4, ...
                num2cell(1:N), ...
                'UniformOutput', false);
        end
        function t4c  = t4FramesCells(~, bname, findices)
            assert(ischar(bname));
            assert(isnumeric(findices));
            t4c = cellfun( ...
                @(x) sprintf('%s_frame%i_to_resolved_t4', bname, x), ...
                num2cell(findices), ...
                'UniformOutput', false);
        end
        function t4s  = t4sInv(this, t4s)
            for t = 1:length(t4s)
                this.buildVisitor.t4_inv(t4s{t});
                t4s{t} = this.buildVisitor.filenameT4inv(t4s{t});
            end
        end
        function t4s  = t4sMul(this, t4s1, t4s2)
            for t = 1:min(length(t4s1), length(t4s2))
                t4s{t} = this.buildVisitor.filenameT4mul(t4s1{t}, t4s2{t});
                this.buildVisitor.t4_mul(t4s1{t}, t4s2{t}, t4s{t});                
            end
        end
        function t4s  = t4sMuls(this, varargin)
            nva = length(varargin);
            t4s = varargin{end};
            for v = nva-1:-1:1
                t4s = this.t4sMul(t4s, varargin{v});
            end
        end
        function        t4sImg(this, t4s, fdfp0, fdfp1, ref_fdfp)
            assert(iscell(t4s));
            assert(iscell(fdfp1) && all(size(t4s) == size(fdfp1)));
            fdfp0    = this.char2cell(fdfp0, size(t4s));
            ref_fdfp = this.char2cell(ref_fdfp, size(t4s));
            for f = 1:length(t4s)
                this.buildVisitor.t4img_4dfp(t4s{f}, fdfp0{f}, 'out', fdfp1{f}, 'options', ['-O' ref_fdfp{f}]);
            end
        end
        function fdfp = zeroSlicesOnBoundaries(this, fdfp, n)
            N     = this.readSize(fdfp, 3);
            fdfp1 = sprintf('%s_z%ito%i', fdfp,  1, n);
            fdfp2 = sprintf('%s_z%ito%i', fdfp1, N-n+1, N);
            this.buildVisitor.zero_slice_4dfp(fdfp,  'z', 1, n,     fdfp1);
            this.buildVisitor.zero_slice_4dfp(fdfp1, 'z', N-n+1, N, fdfp2);
            fdfp  = fdfp2;            
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

 