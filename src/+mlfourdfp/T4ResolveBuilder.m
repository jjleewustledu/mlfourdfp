classdef T4ResolveBuilder < mlfourdfp.AbstractT4ResolveBuilder
	%% T4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 10-Mar-2016 21:29:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.  Copyright 2017 John Joowon Lee.
    
    
	methods
 		function this = T4ResolveBuilder(varargin)
 			%% T4RESOLVEBUILDER
            %  @param theImages =: safe fileprefix =: time summed.
            %  @param blurArg; default := this.sessiondata.compositeT4ResolveBuilderBlurArg.
            %  @param indicesLogical 
            %  @param indexOfReference; default := 1, the first of theImages.
 			
            this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addOptional( ip, 'cctor', []);
            addParameter(ip, 'blurArg', this.sessionData.t4ResolveBuilderBlurArg, @isnumeric);
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});
            
            if (isempty(ip.Results.cctor))
                import mlfourdfp.*;
                this.imageComposite_ = ImageFrames( ...
                    this, ...
                    'theImages', this.ensureSafeFileprefix(this.theImages), ...
                    'indicesLogical', ip.Results.indicesLogical, ...
                    'indexOfReference', ip.Results.indexOfReference);
                this.blurArg_ = ip.Results.blurArg;
            end
            this = this.updateFinished;
        end
                
        function this         = resolve(this, varargin)
            import mlfourdfp.*;
            ip = inputParser;
            addParameter(ip, 'dest',           '',                  @ischar);
            addParameter(ip, 'source',         this.theImages,      @FourdfpVisitor.lexist_4dfp);
            addParameter(ip, 'destMask',       'none',              @ischar);
            addParameter(ip, 'sourceMask',     'none',              @ischar);
            addParameter(ip, 'destBlur',       this.blurArg,        @isnumeric); % fwhh/mm
            addParameter(ip, 'sourceBlur',     this.blurArg,        @isnumeric); % fwhh/mm 
            addParameter(ip, 'maskForImages',  'Msktgen',           @ischar);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 't40',            this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'resolveTag',     this.resolveTag,     @ischar);
            addParameter(ip, 'log',            '/dev/null',         @ischar);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            this.resolveTag = ip.Results.resolveTag;            
            ipr = ip.Results;
            ipr = this.expandBlurs(ipr);       
            ipr = this.expandMasks(ipr);
            ipr.source = this.ensureLocalFourdfp(ipr.source);
            if (isempty(ipr.dest)); ipr.dest = ipr.source; end
            ipr.resolved = ipr.source; % initialize this.revise   
            if (~any(this.indicesLogical))
                ipr.dest = this.fileprefixRevision(ipr.dest, this.NRevisions);            
                ipr.resolved = sprintf('%s_%s', ipr.dest, this.resolveTag);
                this.buildVisitor_.copyfile_4dfp(ipr.source, ipr.resolved);
                this = this.finalize(ipr);
                return
            end
            if (this.isfinished)  
                this = this.alreadyFinalized(ipr);
                return
            end
            while (this.rnumber <= this.NRevisions)
                ipr.source = ipr.resolved;
                ipr.dest   = this.fileprefixRevision(ipr.dest, this.rnumber);
                [ipr,this] = this.revise(ipr);
                assert(this.rnumber < 10);
            end
            this = this.finalize(ipr);
        end
        function [ipr,this]   = revise(this, ipr)
            this.copySourceToDest(ipr); % crop/copy ipr.source to ipr.dest
            this.imageRegLog = loggerFilename( ...
                ipr.dest, 'func', 'T4ResolveBuilder_imageReg', 'path', this.logPath);
            this.resolveLog = loggerFilename( ...
                ipr.dest, 'func', 'T4ResolveBuilder_resolveAndPaste', 'path', this.logPath);
            
            this = this.imageReg(ipr);
            ipr = this.resolveAndPaste(ipr); 
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end
        function this =         imageReg(this, ipr)
            stagedImgs  = this.lazyStageImages(ipr);                 % contracted wrt this.indicesLogical
            blurredImgs = this.lazyBlurImages(ipr);                  % "
            maskedImgs  = this.lazyMasksForImages(ipr, blurredImgs); % "
            assertSizeEqual(stagedImgs, blurredImgs, maskedImgs);
            len = sum(this.indicesLogical);
            t4Failures = zeros(len, len);
            for m = 1:len
                for n = 1:len
                    if (m ~= n) 
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            if (~lexist(t4))
                                this.buildVisitor.(this.sessionData.frameAlignMethod)( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
                                    'destMask',   maskedImgs{m}, ...
                                    'sourceMask', maskedImgs{n}, ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'log',        this.imageRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            t4Failures(m,n) = t4Failures(m,n) + 1;
                            copyfile( ...
                                this.buildVisitor.transverse_t4, ...
                                this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m}), 'f');
                            dispwarning(ME);
                        end
                    end
                end
            end
            t4Failures = sum(t4Failures, 1);
            fprintf('T4ResolveBuilder.imageReg: t4Failures->%s\n', mat2str(t4Failures));
            this.indicesLogical(this.indicesLogical) = ...
                ensureRowVector(this.indicesLogical(this.indicesLogical)) & ...
                ensureRowVector(t4Failures < 0.25*len);
            fprintf('T4ResolveBuilder.imageReg: this.indicesLogical->%s\n', mat2str(this.indicesLogical)); 
            
            this.t4_resolve_err = nan(len, len);
            for m = 1:length(stagedImgs)
                for n = 1:length(stagedImgs)
                    if (m ~= n)
                        try
                            [rmsdeg,rmsmm] = this.t4_resolve_errParser(this.resolvePair( ...
                                    mybasename(stagedImgs{m}), mybasename(stagedImgs{n})));
                            this.t4_resolve_err(m,n) = this.t4_resolve_errAverage(rmsdeg, rmsmm);
                        catch ME
                            dispwarning(ME);
                        end
                        
                    end
                end
            end 
            fprintf('T4ResolveBuilder.imageReg: this.t4_resolve_err->%s\n', mat2str(this.t4_resolve_err));
            
            this.deleteTrash;
        end
        function [ipr,imgFps] = resolveAndPaste(this, ipr)
            %% RESOLVEANDPASTE - preassign ipr.dest, this.resolveTag, this.indexOfReference as needed.
            %  @param ipr is a struct w/ field dest, a string fileprefix || is a string.

            assert(isstruct(ipr));
            
            pwd0   = pushd(fileparts(ipr.dest));
            imgFps = mybasename(this.fileprefixOfReference(ipr)); % initial on ipr.dest
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    %                    fileprefix of frame != this.indexOfReference
                    ipr.currentIndex = f;
                    imgFps = [imgFps ' ' mybasename(this.fileprefixIndexed(ipr))]; %#ok<AGROW>
                end
            end  
            %% Must use short fileprefixes in calls to t4_resolve to avoid filenaming error by t4_resolve  
            %  t4_resolve: /data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/V1/FDG_V1-NAC/E8/fdgv1e8r1_frame1_to_/data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/V1/FDG_V1-NAC/E8/fdgv1e8r1_frame8_t4 read error          
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, imgFps, ...
                'options', '-v -m -s', 'log', this.resolveLog);
            this.t4imgAll(ipr, this.resolveTag); % transform ipr.dest on this.resolveTag
            this.reconstituteImages(ipr, this.resolveTag); % reconstitute all frames            
            ipr.resolved = sprintf('%s_%s', ipr.dest, this.resolveTag);
            popd(pwd0);
        end
        function                reconstituteImages(this, ipr, varargin)
            if (this.skipT4imgAll)
                return
            end
            
            ip = inputParser;
            addRequired(ip, 'ipr', @isstruct);
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, ipr, varargin{:});            
            tag = mybasename(ip.Results.tag);
            
            import mlpet.*;
            prev = PETImagingContext([ipr.dest '.4dfp.ifh']);
            prev = prev.numericalNiftid; 
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    curr = PETImagingContext([this.fileprefixIndexedResolved(ipr.dest, f, tag) '.4dfp.ifh']);
                    curr = curr.numericalNiftid;
                    prev.img(:,:,:,f) = curr.img(:,:,:);
                end
            end
            if (~isempty(tag))
                prev.filename = [ipr.dest '_' tag '.4dfp.ifh'];
            else
                prev.filename = [ipr.dest '.4dfp.ifh'];
            end
            prev.save;
            indicesLogical = this.indicesLogical; %#ok<NASGU>
            save([prev.fqfileprefix '_indicesLogical.mat'], 'indicesLogical');
        end
        function                teardownResolve(this, ipr)
            try
                for r = 1:this.NRevisions                
                    fp0 = this.fileprefixRevision(ipr.dest, r);
                    for il = 1:length(this.indicesLogical)
                        if (this.indicesLogical(il))
                            deleteExisting(sprintf('%s_frame%i.4dfp.*',    fp0, il));
                            deleteExisting(sprintf('%s_frame%i_b*.4dfp.*', fp0, il));
                            deleteExisting(sprintf('%s_frame%i_C*.4dfp.*', fp0, il));
                            deleteExisting(sprintf('%s_frame%i_%s.4dfp.*', fp0, il, this.resolveTag));
                            deleteExisting(sprintf('%s_frame%i_g*.4dfp.*', fp0, il));
                            deleteExisting(sprintf('%s_frame%i_g*.nii.gz', fp0, il));
                        end
                    end
                end            
                deleteExisting(sprintf('%s_*_*.*', ipr.maskForImages));
                deleteExisting(sprintf('*_g0_1.4dfp.*'));
            catch ME
                handwarning(ME);
            end
        end
        
        %% UTILITY
              
        function         copySourceToDest(this, ipr)
            if (this.skipT4imgAll)
                return
            end
            
            if (1 == this.rnumber)
                try
                    if (~this.buildVisitor.lexist_4dfp(ipr.dest) && ...
                        ~strcmp(ipr.source, ipr.dest))
                        this.buildVisitor.copy_4dfp(ipr.source, ipr.dest);
                    end
                catch ME
                    handexcept(ME);
                end
                return
            end
            this.copyResolvedToNewRevision(ipr);
        end
        function         copyResolvedToNewRevision(this, ipr)
            %% COPYRESOLVEDTONEWREVISION opportunistically reuses existing files from the last iteration
            
            this.buildVisitor.copy_4dfp( ...
                this.fileprefixResolved(ipr.dest, this.rnumber-1), ...
                this.fileprefixRevision(ipr.dest, this.rnumber));
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
        function fp    = fileprefixIndexedResolved(this, varargin)
            ip = inputParser;
            addRequired(ip, 'dest', @ischar);
            addRequired(ip, 'currentIndex', @isnumeric);
            addOptional(ip, 'tag', this.resolveTag, @ischar);
            parse(ip, varargin{:});
            ipr = ip.Results;
            fp = sprintf('%s_%s', this.fileprefixIndexed(ipr), ip.Results.tag);
        end
        function fp    = fileprefixOfReference(this, ipr, varargin)
            ip = inputParser;
            addOptional(ip, 'indexOfRef', this.indexOfReference, @isnumeric)
            parse(ip, varargin{:});            
            ipr.currentIndex = ip.Results.indexOfRef;            
            fp = this.fileprefixIndexed(ipr); 
        end
        function fqfps = lazyMasksForImages(this, ipr, blurredImgs)
            %  @param ipr.maskForImages is the base fileprefix for masks.
            %  @param ipr.sourceMask is cell-array of fileprefixes for masks or cells of 'none'.
            %  @param ipr.destmask   is cell-array of fileprefixes for masks or cells of 'none'.
            %  @param blurredImgs are f.-q.-fileprefixes images to mask.
            %  @return fqfps is a cell-array of f.-q. mask fileprefixes.
            
            fqfps = cell(1, length(blurredImgs));            
            if (isempty(ipr.maskForImages) || strcmpi(ipr.maskForImages, 'none'))
                fqfps = cellfun(@(x) 'none', fqfps, 'UniformOutput', false); 
                return
            end            
              
            if (strcmp(ipr.maskForImages, 'Msktgen'))
                mg    = mlpet.Msktgen('sessionData', this.sessionData);
                mskt  = mg.constructMskt( ...
                    'source', this.sessionData.tracerRevision, ...
                    'intermediaryForMask', this.sessionData.T1001, ...
                    'sourceOfMask', fullfile(this.sessionData.vLocation, 'brainmask.4dfp.ifh'), ...
                    'blurForMask', 6*this.blurArg, 'threshp', 0);
                mskt.save;
                fqfps = cellfun(@(x) mskt.fqfileprefix, fqfps, 'UniformOutput', false);
                return
            end
            
            % build masks for each frame
            sourceSumt = this.buildSourceTimeSummed(ipr);
            for b = 1:length(blurredImgs)
                fqfps{b} = sprintf('%s_%i', ipr.maskForImages, b);
                if (~lexist_4dfp(fqfps{b}))
                    this.buildMaskAdjustedForImage(fqfps{b}, blurredImgs{b}, ipr.sourceMask{b}, sourceSumt);
                end
            end
        end
        function fqfps = lazyStageImages(this, ipr)
            fqfps = this.imageComposite.lazyExtractImages(ipr);
        end
        function fqfp  = buildSourceTimeSummed(~, ipr)
            %  @param ipr.source is a f.-q.-fileprefix.
            %  @return fqfp := [ipr.source '_sumt'] generated on the filesystem.  
            %  See also mlfourd.ImagingContext.timeSummed.            
            if (lexist_4dfp([ipr.source '_sumt']))
                fqfp = [ipr.source '_sumt'];
                return
            end
            
            ic = mlfourd.ImagingContext(ipr.source);
            ic = ic.timeSummed;
            ic.save;
            fqfp = ic.fqfileprefix;
        end
        function         buildMaskAdjustedForImage(this, fqfp, blurredImg, sourceMask, sourceSumt)
            bv = this.buildVisitor;
            if (strcmpi(sourceMask, 'msktgen_4dfp'))
                bv.msktgen_4dfp(sourceSumt, 0);
                bv.move_4dfp([sourceSumt '_mskt'], fqfp);
                return
            end
            if (strcmpi(sourceMask, 'msktgen_b110_4dfp'))
                bv.msktgen_b110_4dfp(sourceSumt, 0);
                bv.move_4dfp([sourceSumt '_mskt'], fqfp);
                return
            end
            if (strcmpi(sourceMask, 'msktgen2_4dfp'))
                bv.msktgen2_4dfp(sourceSumt, 0);
                bv.move_4dfp([sourceSumt '_mskt'], fqfp);
                return
            end  
            this.buildMaskForImage(blurredImg, fqfp); % more efficient than calling buildMasksForImages
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function this = t4imgAll(this, ipr, tag)
            if (this.skipT4imgAll) 
                return
            end
            tag = mybasename(tag);
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    frameFp = sprintf('%s_frame%i', ipr.dest, f);
                    this.buildVisitor.t4img_4dfp( ...
                        sprintf('%s_to_%s_t4', frameFp, tag), ...
                        frameFp, ...
                        'out', sprintf('%s_%s', frameFp, tag), ...
                        'options', ['-O' frameFp]);
                end
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
    
    %% HIDDEN @deprecated
    
    methods (Hidden)
        function [mskm,mskn] = lazyMaskForImages_deprecated(this, ipr, fpm, fpn, m, n)
            if (isempty(ipr.maskForImages))
                mskm = 'none';
                mskn = 'none';
                return
            end
            
            maskMN = sprintf('%s_%i_%i.4dfp.img', ipr.maskForImages, m, n);
            maskNM = sprintf('%s_%i_%i.4dfp.img', ipr.maskForImages, n, m);
            if (lexist(maskMN, 'file'))
                mskm = maskMN;
                mskn = mskm;
            else
                if (lexist(maskNM, 'file'))
                    mskm = maskNM;
                    mskn = mskm;
                else
                    mskm = this.buildMaskForImages(fpm, fpn, maskMN);
                    mskn = mskm;                    
                end
            end
        end
    end
end

 