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
 			
            this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addOptional( ip, 'cctor', []);
            addParameter(ip, 'blurArg', this.sessionData.t4ResolveBuilderBlurArg, @isnumeric);
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'theImages', {}, @(x) iscell(x) || ischar(x));
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});
            ims = ensureCell(ip.Results.theImages);
            assert(lexist(ims{1}) || this.buildVisitor.lexist_4dfp(ims{1}), ....
                'mlfourdfp.T4ResolveBuilder.ctor.ip.Results.theImages->%s does not exist', ...
                cell2str(ims));
            
            if (isempty(ip.Results.cctor))
                import mlfourdfp.*;
                this.imageComposite_ = ImageFrames(this, ...
                    'indicesLogical', ip.Results.indicesLogical, ...
                    'theImages', FourdfpVisitor.ensureSafeFileprefix(ip.Results.theImages), ...
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
            addParameter(ip, 'maskForImages',  'maskForImages',     @ischar);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 't40',            this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'resolveTag',     this.resolveTag,     @ischar);
            addParameter(ip, 'log',            '/dev/null',         @ischar);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            this.resolveTag = ip.Results.resolveTag;            
            ipr = ip.Results;
            ipr = this.expandBlurs(ipr);            
            if (isempty(ipr.dest)); ipr.dest = ipr.source; end
            ipr.resolved = ipr.source; % initialize this.revise   
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
            
            this.imageReg(ipr);
            ipr = this.resolveAndPaste(ipr); 
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end
        function                imageReg(this, ipr)
            stagedImgs  = this.lazyStageImages(ipr);
            blurredImgs = this.lazyBlurImages(ipr);
            assert(length(stagedImgs) == length(blurredImgs));
            for m = 1:length(stagedImgs)
                for n = 1:length(stagedImgs)
                    if (m ~= n) 
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            if (~lexist(t4))
                                maskFp = this.lazyMaskForImages( ...
                                    ipr.maskForImages, stagedImgs{m}, stagedImgs{n}, ...
                                    this.imageComposite.fortranImageIndices(m), this.imageComposite.fortranImageIndices(n));
                                this.buildVisitor.(this.sessionData.alignMethod)( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
                                    'destMask',   maskFp, ...
                                    'sourceMask', maskFp, ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'log',        this.imageRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            copyfile( ...
                                this.buildVisitor.transverse_t4, ...
                                this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m}), 'f');
                            handwarning(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end
        function [ipr,imgFps] = resolveAndPaste(this, ipr)
            %% RESOLVEANDPASTE - preassign ipr.dest, this.resolveTag, this.indexOfReference as needed.
            %  @param ipr is a struct w/ field dest, a string fileprefix || is a string.

            assert(isstruct(ipr));
            
            pwd0   = pushd(fileparts(ipr.dest));
            imgFps = mybasename(this.fileprefixOfReference(ipr, this.indexOfReference)); % initial on ipr.dest
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    %                    fileprefix of frame != this.indexOfReference
                    imgFps = [imgFps ' ' mybasename(this.fileprefixIndexed(ipr.dest, f))]; %#ok<AGROW>
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
        end             
        function this         = finalize(this, ipr)
            this.ipResults_ = ipr;
            this.rnumber = this.NRevisions;
            this.product_ = mlpet.PETImagingContext([ipr.resolved '.4dfp.ifh']);            
            %assert(~isempty(this.product_));
            this.teardownResolve(ipr);
            this.finished.touchFinishedMarker;          
        end
        function this         = alreadyFinalized(this, ipr)
            dest = this.fileprefixRevision(ipr.dest, this.NRevisions);
            ipr.resolved = sprintf('%s_%s', dest, this.resolveTag);            
            this.ipResults_ = ipr;
            this.rnumber = this.NRevisions;
            this.product_ = mlpet.PETImagingContext([ipr.resolved '.4dfp.ifh']);
            %assert(~isempty(this.product_));
        end
        function                teardownResolve(this, ipr)
            %if (this.keepForensics); return; end
            
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
        function fp    = fileprefixIndexed(~, fp, fr)
            assert(ischar(fp));
            assert(isnumeric(fr));
            fp = sprintf('%s_frame%i', fp, fr);
        end
        function fp    = fileprefixIndexedResolved(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fp', @ischar);
            addRequired(ip, 'fr', @isnumeric);
            addOptional(ip, 'tag', this.resolveTag, @ischar);
            parse(ip, varargin{:});
            fp = sprintf('%s_%s', this.fileprefixIndexed(ip.Results.fp, ip.Results.fr), ip.Results.tag);
        end
        function fp    = fileprefixOfReference(this, ipr, varargin)
            ip = inputParser;
            addOptional(ip, 'indexOfRef', this.indexOfReference, @isnumeric)
            parse(ip, varargin{:});
            
            fp = this.fileprefixIndexed(ipr.dest, ip.Results.indexOfRef); % mybasename =: BUG?
        end
        function fqfp  = lazyBlurImage(this, ipr, blur)
            %% LAZYBLURIMAGE uses specifiers in ipr; will not replace any existing image
            %  @param ipr is a struct
            %  @param blur is numeric
            
            fqfp_ = this.fileprefixIndexed(ipr.dest, ipr.currentIndex);
            fqfp  = this.fileprefixBlurred(fqfp_, blur);
            if (~this.buildVisitor.lexist_4dfp(fqfp))
                this.buildVisitor.imgblur_4dfp(fqfp_, blur);
            end
        end
        function fp    = lazyMaskForImages(this, maskFp, fpm, fpn, m, n)
            if (isempty(maskFp))
                fp = '';
                return
            end
            
            maskMN = sprintf('%s_%i_%i.4dfp.img', maskFp, m, n);
            maskNM = sprintf('%s_%i_%i.4dfp.img', maskFp, n, m);
            if (lexist(maskMN, 'file'))
                fp = maskMN;
            else
                if (lexist(maskNM, 'file'))
                    fp = maskNM;
                else
                    fp = this.maskForImages(fpm, fpn, maskMN);
                end
            end
        end
        function fqfps = lazyStageImages(this, ipr)
            fqfps = this.imageComposite.lazyExtractImages(ipr);
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function this = t4imgAll(this, ipr, tag)
            if (this.skipT4imgAll) % || this.rnumber < this.NRevisions)
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
end

 