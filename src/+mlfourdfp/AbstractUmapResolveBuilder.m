classdef (Abstract) AbstractUmapResolveBuilder < mlfourdfp.CompositeT4ResolveBuilder
	%% ABSTRACTUMAPRESOLVEBUILDER  

	%  $Revision$
 	%  was created 14-Nov-2016 22:37:22
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	
    
    properties (Constant)
        REPLACE_COMPLETED = true
    end

	properties
        ct_rescaleSlope = 1
        ct_rescaleIntercept = -1024
        reuseCTMasked = true
        reuseCTRescaled = true
        reuseCarneyUmap = true
        sessionDataCache
    end
    
	methods
        
 		function this  = AbstractUmapResolveBuilder(varargin)
 			%% ABSTRACTUMAPRESOLVEBUILDER
 			%  Usage:  this = AbstractUmapResolveBuilder()

 			this = this@mlfourdfp.CompositeT4ResolveBuilder(varargin{:});
            this.mmrBuilder_ = mlsiemens.MMRBuilder('sessionData', this.sessionData);
            this.NRevisions = 2;
            this.blurArg_ = 1.5; % per Avi, 2016oct25
        end
        
        function this  = buildUmap(this, varargin)
        end
        function this  = buildUmaps(this, varargin)
        end        
        function dest  = buildTracerNAC(this, varargin)
            %% BUILDTRACERNAC builds 4dfp formatted NAC images.
            %  See also:  mlfourdfp.FourdfpVisitor.sif_4dfp.
            
            ip = inputParser;
            addOptional(ip, 'sessionData', this.sessionData, @(x) isa(x, 'mlpipeline.ISessionData'));
            parse(ip, varargin{:});            
            sessd = ip.Results.sessionData;
            sessd.attenuationCorrected = false;
            lm      = sessd.tracerListmodeSif(      'typ', 'fqfp');
            dest    = sessd.tracerRevision('typ', 'fqfp');
            destLoc = sessd.tracerRevision('typ', 'path');
            
            if (this.buildVisitor.lexist_4dfp(dest))
                return
            end
            if (~this.buildVisitor.lexist_4dfp(lm))
                fprintf('mlfourdfp.AbstractUmapResolveBuilder.buildTracerNAC.buildVisitor.sif_4dfp is building %s\n', lm);
                pwd0 = pwd;
                cd(fileparts(lm));
                this.buildVisitor.sif_4dfp(lm);
                cd(pwd0);
            end
            if (~isdir(destLoc))
                mkdir(destLoc);
            end
            ipr = struct('dest', dest, 'source', lm, 'rnumber', 0);
            this.copySourceToDest(ipr);
        end
        function this  = convertUmapTo4dfp(this)
            pwd0 = pwd;
            if (~lexist(this.sessionData.tracerListmodeUmap('typ', '4dfp.img')))
                cd( this.sessionData.tracerListmodeUmap('typ', 'path'));
                this.buildVisitor.IFhdr_to_4dfp( ...
                    this.sessionData.tracerListmodeUmap('typ', 'v.hdr'));
                cd(pwd0);
            end
        end
        function this  = convertUmapToE7Format(this, varargin)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            
            ip = inputParser;
            addOptional(ip, 'umap', ...
                fullfile(sessd.tracerNACLocation, ...
                    sprintf('%s_op_%s', sessd.umapSynth('typ', 'fp'), sessd.tracerNACRevision('typ', 'fp'))), ...
                @lexist_4dfp);
            addParameter(ip, 'zoom', this.mmrBuilder_.inverseCrop, @isnumeric);
            parse(ip, varargin{:});
            umap = ip.Results.umap;
            
            flipped = this.buildVisitor.flip_4dfp('z', umap);
            ic = mlfourd.ImagingContext([flipped '.4dfp.ifh']);
            ic = ic.zoomed(ip.Results.zoom);
            ic.noclobber = false;
            ic.saveas([flipped '.4dfp.ifh']);
            movefile( ...
                sprintf('%s.4dfp.img', flipped), ...
                sprintf('%s.v',        umap), 'f');
            if (~this.keepForensics)
                delete(sprintf('%s.4dfp.*', flipped));
                delete(sprintf('%sfz.4dfp.*', umap));
                delete(sprintf('%s*.log', umap));
            end
            
            this.product_ = mlfourd.ImagingContext(sprintf('%s.v', umap));
        end
        function this  = convertUmapsToE7Format(this, umaps)
            assert(iscell(umaps));
            prodCell = {};
            for fr = 1:length(umaps)
                this = this.convertUmapToE7Format(umaps{fr});
                prodCell{fr} = this.product_.fqfilename; %#ok<AGROW>
            end
            
            this.product_ = prodCell;
        end  
        function this  = repUmapToE7Format(this, umaps)
            assert(iscell(umaps));
            prodCell = {};
            for fr = 1:length(umaps)
                this = this.repUmapToE7Format__(umaps{fr});
                prodCell{fr} = this.product_.fqfilename; %#ok<AGROW>
            end
            
            this.product_ = prodCell;
        end 
        function this  = repUmapToE7Format__(this, varargin)
            sessd = this.sessionData;
            sessd.rnumber = 1;
            
            ip = inputParser;
            addOptional(ip, 'umap', ...
                fullfile(sessd.tracerNACLocation, ...
                    sprintf('%s_op_%s', sessd.umapSynth('typ', 'fp'), sessd.tracerNACRevision('typ', 'fp'))), ...
                @lexist_4dfp);
            addParameter(ip, 'zoom', this.mmrBuilder_.inverseCrop, @isnumeric);
            parse(ip, varargin{:});
            umap = ip.Results.umap;
            
            flipped = this.buildVisitor.flip_4dfp('z', umap);
            ic = mlfourd.ImagingContext([flipped '.4dfp.ifh']);
            ic = ic.zoomed(ip.Results.zoom);
            ic.noclobber = false;
            ic.saveas([flipped '.4dfp.ifh']);
            movefile( ...
                sprintf('%s.4dfp.img', flipped), ...
                sprintf('%s.v',        umap), 'f');
            if (~this.keepForensics)
                delete(sprintf('%s.4dfp.*', flipped));
                delete(sprintf('%sfz.4dfp.*', umap));
                delete(sprintf('%s*.log', umap));
            end
            
            this.product_ = mlfourd.ImagingContext(sprintf('%s.v', umap));
        end 
        function this  = loadSessionDataCache(this, varargin)
            ip = inputParser;
            addOptional(ip, 'tracers', {'OC' 'HO' 'OO'}, @iscell);
            parse(ip, varargin{:});
            
            cache = {};
            for s = 1:2
                for t = 1:length(ip.Results.tracers)
                    sessd = this.sessionData;
                    sessd.tracer = ip.Results.tracers{t};
                    sessd.snumber = s;
                    if (isdir(sessd.tracerListmodeLocation))
                        cache = [cache {sessd}]; %#ok<AGROW>
                    end
                end
            end
            this.sessionDataCache = cache;
        end
        function         reconvertUmapsToE7Format(this)
            pwd0 = pushd(this.sessionData.fdgNACLocation);
            for fr = 1:length(this.indicesLogical)
                this.buildVisitor.extract_frame_4dfp(this.umapsOpTracer, fr);
                this.convertUmapToE7Format(this.umapsOpTracerFrame(fr));
                delete([this.umapsOpTracerFrame(fr) '.4dfp.*']);
                delete([this.umapsOpTracerFrame(fr) '_flipz.log']);
            end
            popd(pwd0);
        end
        function loc   = resolveSequenceLocation(this, varargin)
            %  @param named tracer is a string identifier.
            %  @param named snumber is the scan number; is numeric.
            %  @param named typ is string identifier:  folder path, fn, fqfn, ...  
            %  See also:  imagingType.
            %  @param named frame is numeric.
            %  @param named rnumber is the revision number; is numeric.
            %  @returns ipr, the struct ip.Results obtained by parse.            
            %  @returns schr, the s-number as a string.
            
            ipr = this.sessionData.iprLocation(varargin{:});
            tag = this.resolveSequenceTag;
            loc = locationType(ipr.typ, ...
                    fullfile(this.sessionData.tracerNACLocation( ...
                        'tracer', this.sessionData.tracer, ...
                        'snumber', this.sessionData.snumber), ...
                    sprintf('%s%s', upper(tag(1)), tag(2:end)), '')); 
        end
        function obj   = resolveSequenceResolved(this, varargin)
            %  @param named tracer is a string identifier.
            %  @param named snumber is the scan number; is numeric.
            %  @param named typ is string identifier:  folder path, fn, fqfn, ...  
            %  See also:  imagingType.
            %  @param named frame is numeric.
            %  @param named rnumber is the revision number; is numeric.
            %  @returns ipr, the struct ip.Results obtained by parse.            
            %  @returns schr, the s-number as a string.
            
            fqfn = fullfile( ...
                this.resolveSequenceLocation, ...
                sprintf('%sr%i_resolved.4dfp.ifh', this.resolveSequenceTag, this.NRevisions));
            obj  = this.sessionData.fqfileprefixObject(fqfn, varargin{:});
        end
        function fp    = resolveSequenceTag(this)
            sessd = this.sessionData;
            if (isempty(sessd.tracer))
                fp = sprintf('umapResolveSequencev%i', sessd.vnumber);
                return
            end
            fp = sprintf('%sUmapResolveSequencev%i', sessd.tracer, sessd.vnumber);
        end
        function umaps = umapsOpTracer(this)
            umaps = sprintf('umapsOp%sv%ir%i', ...
                upperFirst(this.sessionData.tracer), this.sessionData.vnumber, this.sessionData.rnumber);
        end
        function fp    = umapsOpTracerFrame(this, fr)
            ipr = struct( ...
                'dest', this.umapsOpTracer, ...
                'currentIndex', fr);
            fp = this.fileprefixIndexed(ipr);
        end
        function [s,r] = viewUmaps(this)
            [s,r] = mlbash(sprintf( ...
                'fslview %s.4dfp.hdr %s.4dfp.hdr', ...
                this.sessionData.tracerNACRevision, ...
                this.umapsOpTracer));
        end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        f18UmapResolveBuilder_
        mmrBuilder_
    end
    
    methods (Access = protected)
        function [this,ct]         = alignCTToSumtResolved(this, varargin)
            ip = inputParser;
            addRequired(ip, 'ct',           @(x) lexist(this.fourdfpImg(x), 'file'));
            addOptional(ip, 'sumtResolved', this.sumTimes(this.sessionData.fdgNACResolved('typ', 'fqfp')), ...
                                            @(x) lexist(this.fourdfpImg(x), 'file'));
            parse(ip, varargin{:});
                  
            this.resolveTag = ['op_' mybasename(ip.Results.sumtResolved)];
            this = this.resolveSequence(ip.Results.sumtResolved, this.sessionData.mpr('typ', 'fqfp'), ip.Results.ct);
            ct   = this.product{3}.fqfp;
        end
        function [this,ctOnMpr1]   = alignCTToMpr(this, varargin)
            ip = inputParser;
            addRequired(ip, 'ct',  @(x) lexist(this.fourdfpImg(x), 'file'));
            addOptional(ip, 'mpr', this.sessionData.mpr('typ', 'fqfp'), ...
                                   @(x) lexist(this.fourdfpImg(x), 'file'));
            parse(ip, varargin{:});
            
            ctOnMpr = this.CT2mpr_4dfp(ip.Results.ct);
            ctOnMpr1 = strrep_on_(ctOnMpr);
            this.buildVisitor.lns_4dfp(ctOnMpr, mybasename(ctOnMpr1));
            this.resolveTag = ['op_' mybasename(ip.Results.mpr)];
            this = this.resolveSequence(ip.Results.mpr, mybasename(ctOnMpr1));
            ctOnMpr1 = this.product{2}.fqfp;
        end
        function [ctm,ic]          = buildCTMasked(this)
            %% CTMASKED
            %  @return ctm := this.sessionData.ctMasked('typ', 'fqfp')
            %  @return ic  := ctMasked as ImagingContext on CT-space
            
            import mlfourd.*;
            mpr  = this.sessionData.mpr('typ', 'fqfp');
            ct   = this.sessionData.ct('typ', 'fqfp');
            ctm  = this.sessionData.ctMasked('typ', 'fqfp');
            if (lexist(this.fourdfpImg(ctm)) && this.reuseCTMasked)
                ic = ImagingContext(ctm);
                return
            end
            
            [ctOnMpr,ctToMprT4] = this.CT2mpr_4dfp(ct, ...
                'log', sprintf('CarneyUmapBuilder_CT2mpr_4dfp_%s.log', datestr(now,30)));
            mprToCtT4 = this.buildVisitor.t4_inv(ctToMprT4);
            mprb = this.buildVisitor.imgblur_4dfp(mpr, 10);
            
            ct_  = sprintf('%s_%s', ct, datestr(now, 30));
            ct_  = this.buildVisitor.maskimg_4dfp(ctOnMpr, mprb, ct_, 'options', '-t5'); % in mpr-space
            ct__ = sprintf('%s%s', ct_, 'a');
            ct__ = this.buildVisitor.maskimg_4dfp(ct_, ct_, ct__, 'options', '-t50');
            ctm  = this.buildVisitor.t4img_4dfp(mprToCtT4, ct__, 'out', ctm, 'options', ['-O' ct]); % back to ct-space
            ic = ImagingContext(ctm);
            delete([ct_ '.4dfp.*']);
            %delete([ct__ '.4dfp.*']); ct__ in mpr-space has best
            %registration
        end
        function [ctm,ic,ctToMprT4] = buildCTMasked2(this)
            %% CTMASKED
            %  @return ctm := this.sessionData.ctMasked('typ', 'fqfp')
            %  @return ic  := ctMasked as ImagingContext on MPR-space
            
            import mlfourd.*;
            mpr  = this.sessionData.mpr('typ', 'fqfp');
            ct   = this.sessionData.ct('typ', 'fqfp');
            ctm  = this.sessionData.ctMasked('typ', 'fqfp');
            if (lexist(this.fourdfpImg(ctm)) && this.reuseCTMasked)
                ic = ImagingContext(ctm);
                return
            end
            
            [ctOnMpr,ctToMprT4] = this.CT2mpr_4dfp(ct, ...
                'log', sprintf('CarneyUmapBuilder_CT2mpr_4dfp_%s.log', datestr(now,30)));
            mprb = this.buildVisitor.imgblur_4dfp(mpr, 10);
            
            ct_  = sprintf('%s_%s', ct, datestr(now, 30));
            this.buildVisitor.maskimg_4dfp(ctOnMpr, mprb, ct_, 'options', '-t5'); % in mpr-space
            this.buildVisitor.maskimg_4dfp(ct_, ct_, ctm, 'options', '-t50');
            ic = ImagingContext(ctm);
            delete([ct_ '.4dfp.*']);
            %delete([ct__ '.4dfp.*']); ct__ in mpr-space has best
            %registration
        end
        function imageOnSumt       = ctOnPetSumt(this, ct, petSumt)
            assert(lexist(this.fourdfpImg(ct)));
            assert(lexist(this.fourdfpImg(petSumt)));
            
            [~,ctToMprT4]     = this.CT2mpr_4dfp(ct);
            [this,petToMprT4] = this.petSumt2mpr(petSumt);
            mprToPetT4        = this.buildVisitor.t4_inv(petToMprT4);
            ctToPetT4         = this.buildVisitor.t4_mul(ctToMprT4, mprToPetT4);
            imageOnSumt       = this.buildVisitor.t4img_4dfp(ctToPetT4, ct, 'options', ['-O' petSumt]);
        end 
        function [ctOnMpr,ctToMprT4] = CT2mpr_4dfp(this, ct, varargin)
            assert(lexist(this.fourdfpImg(ct), 'file')); % not necessarily in pwd            
            mpr       = this.sessionData.mpr('typ', 'fqfp');
            pth       = fileparts(mpr);
            ctToMprT4 = fullfile(pth, this.buildVisitor.filenameT4(mybasename(ct), mybasename(mpr))); 
            ctOnMpr   = fullfile(pth, [mybasename(ct) '_on_' mybasename(mpr)]);
            
            if (~lexist(this.fourdfpImg(ctOnMpr)))
                ctOnMpr = this.buildVisitor.CT2mpr_4dfp(mpr, ct, 'options', ['-T' this.atlas('typ', 'fqfp')], varargin{:});
            end
            assert(lexist(ctToMprT4, 'file'));        
        end 
        function imageOnSumt       = mprOnPetSumt(this, mpr, petSumt)
            assert(lexist(this.fourdfpImg(mpr)));
            assert(lexist(this.fourdfpImg(petSumt)));
            
            [this,petToMprT4] = this.petSumt2mpr(petSumt);
            mprToPetT4        = this.buildVisitor.t4_inv(petToMprT4);
            imageOnSumt       = this.buildVisitor.t4img_4dfp(mprToPetT4, mpr, 'options', ['-O' petSumt]);
        end   
        function [this,petToMprT4] = petSumt2mpr(this, varargin)
            ip = inputParser;
            addRequired( ip, 'pet', @lexist_4dfp);
            addParameter(ip, 'pharynx', false, @islogical);
            parse(ip, varargin{:});
            pet = ip.Results.pet;
            
            mpr        = this.sessionData.mpr('typ', 'fqfp');
            petToMprT4 = [pet '_to_' mybasename(mpr) '_t4'];
            if (lexist(petToMprT4, 'file'))
                delete(petToMprT4);
            end
            
            this.buildVisitor.t4_inv( ...
                this.buildVisitor.sagittal_t4, 'out', petToMprT4);
            
            petToMprT4 = this.buildVisitor.align_multiSpectral( ...
                'dest', mpr, ...
                'source', pet, ...
                'destMask', 'none', ...
                'sourceMask', 'none', ...
                'destBlur', this.blurArg, ...
                'sourceBlur', this.blurArg, ...
                't40', petToMprT4, ...
                't4', petToMprT4);              
            mprToAtlT4 = [mybasename(mpr) '_to_' this.atlas('typ', 'fp') '_t4'];
            if (~lexist(mprToAtlT4, 'file'))     
                this = this.buildVisitor.msktgenMprage(mpr, this.atlas('typ', 'fp'));          
            end
            
            petToAtlT4 = [pet '_to_' this.atlas('typ', 'fp') '_t4'];
            this.buildVisitor.t4_mul(petToMprT4, mprToAtlT4, petToAtlT4);            
            this.buildVisitor.t4img_4dfp( ...
                petToAtlT4, pet, ...
                'out', [mybasename(pet), '_on_' this.atlas('typ', 'fp')], ...
                'options', ['-O' this.atlas('typ', 'fqfp')]);            
            this.buildVisitor.msktgen2_4dfp( ...
                pet, this.msktgenThresh, ...
                'options', ['-T' this.atlas('typ', 'fqfp')]);            
        end    
        function ctOut             = rescaleCT(this, varargin)
            ip = inputParser;
            addRequired( ip, 'ctMasked', @lexist_4dfp);
            addParameter(ip, 'ctOut', this.sessionData.ctRescaled('typ', 'fqfp'), @ischar);
            parse(ip, varargin{:});
            ctOut = ip.Results.ctOut;
            if (lexist([ctOut '.4dfp.ifh'], 'file') && this.reuseCTRescaled)
                return
            end
            
            ic = mlfourd.ImagingContext([ip.Results.ctMasked '.4dfp.ifh']);
            ic = ic.numericalNiftid;
            ic = ic * this.ct_rescaleSlope + this.ct_rescaleIntercept;            
            ic.noclobber = false;
            ic.saveas([ctOut '.4dfp.ifh']);
        end  
        function this              = resolveSequence(this, varargin)
            import mlfourdfp.*;
            basenames = cellfun(@(x) mybasename(x), varargin, 'UniformOutput', false);
            for b = 1:length(basenames)
                if (~FourdfpVisitor.lexist_4dfp(basenames{b}))
                    FourdfpVisitor.lns_4dfp(varargin{b});
                end
            end
            this.imageComposite = ImageComposite( ...
                this.imageComposite.it4ResolveBuilder, ...
                'theImages', basenames, ...
                'indicesLogical', true(size(basenames)));
            this = this.resolve( ...
                'source', basenames, ...
                'indicesLogical', true(1, length(basenames)));
        end
        function tr                = tracerForNAC(this)
            tr = sprintf('%sForNAC', lower(this.sessionData.tracer));
        end 
        function fp                = tracerNACRevision_fp_r(this, r)
            fp = sprintf('%sv%ir%i', lower(this.sessionData.tracer), this.sessionData.vnumber, r);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

