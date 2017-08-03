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
            addParameter(ip, 'blurArg', 5.5, @isnumeric);
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'theImages', {}, @(x) iscell(x) || ischar(x));
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});
            
            import mlfourdfp.*;
            this.imageComposite_ = ImageFrames(this, ...
                'indicesLogical', ip.Results.indicesLogical, ...
                'theImages', FourdfpVisitor.ensureSafeFileprefix(ip.Results.theImages), ...
                'indexOfReference', ip.Results.indexOfReference);
            this.blurArg_ = ip.Results.blurArg;
            this.finished = mlpipeline.Finished(this, ...
                'path', this.logPath, 'tag', lower(this.sessionData.tracerRevision('typ','fp')));
            cd(this.sessionData.tracerLocation);
        end
                
        function this         = resolve(this, varargin)
            %% RESOLVE iteratively calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param destMask   "
            %  @param source     "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param maskForImages is a f.q. fileprefix.
            %  @param indicesLogical is logical.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param resolveTag is char.
            %  @param log        is the f.q. filename of the log file.
            
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
                return
            end
            while (this.sessionData.rnumber <= this.NRevisions)
                ipr.source = ipr.resolved;
                ipr.dest   = this.fileprefixRevision(ipr.dest, this.sessionData.rnumber);
                [ipr,this] = this.revise(ipr);
                assert(this.sessionData.rnumber < 10);
            end
            this = this.finalize(ipr);
        end
        function [ipr,this]   = revise(this, ipr)
            this.copySourceToDest(ipr); % crop/copy ipr.source to ipr.dest
            this.imageRegLog = loggerFilename( ...
                ipr.dest, 'func', 'T4ResolveBuilder_imageReg', 'path', this.logPath);
            this.resolveLog = loggerFilename( ...
                ipr.dest, 'func', 'T4ResolveBuilder_t4ResolveAndPaste', 'path', this.logPath);
            
            this.imageReg(ipr);
            ipr = this.resolveAndPaste(ipr); 
            this.teardownRevision;
            this.sessionData_.rnumber = this.sessionData_.rnumber + 1;
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
                                this.buildVisitor.align_2051( ...
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
            %% RESOLVEANDPASTE
            %  @param ipr is a struct w/ field dest, a string fileprefix || is a string 
            
            if (ischar(ipr))
                ipr = struct('dest', ipr);
            end
            dest_ = mybasename(ipr.dest);
            imgFps = this.fileprefixOfReference(ipr);
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    imgFps = [imgFps ' ' this.fileprefixIndexed(dest_, f)]; %#ok<AGROW>
                end
            end                     
            
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, imgFps, ...
                'options', '-v -m -s', 'log', this.resolveLog);
            this.t4imgAll(ipr, this.resolveTag);
            this.reconstituteImages(ipr, this.resolveTag);
            
            ipr.resolved = sprintf('%s_%s', dest_, this.resolveTag); 
            %movefile([this.resolveTag '.mat0'], [ipr.resolved '_' datestr(now, 30) '.mat0']);
            %movefile([this.resolveTag '.sub'],  [ipr.resolved '_' datestr(now, 30) '.sub']);
        end
        function                reconstituteImages(this, ipr, varargin)
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
            this.buildVisitor.imgblur_4dfp(ipr.resolved, this.blurArg);
            this.product_ = mlpet.PETImagingContext([ipr.resolved '.4dfp.ifh']);
            this.teardownResolve(ipr);
            this.finished.touchFinishedMarker;            
        end
        function                teardownResolve(this, ipr)
            if (this.keepForensics); return; end
            
            for r = 1:this.NRevisions                
                fp0 = this.fileprefixRevision(ipr.dest, r);
                for f = 1:length(this.indicesLogical)
                    if (this.indicesLogical(f))
                        %delete(sprintf('%s_frame%i.4dfp.*', fp0, f));
                        delete(sprintf('%s_frame%i_b*.4dfp.*', fp0, f));
                        delete(sprintf('%s_frame%i_C*.4dfp.*', fp0, f));
                        %delete(sprintf('%s_frame%i_%s.4dfp.*', fp0, f, this.resolveTag));
                        delete(sprintf('%s_frame%i_g*.nii.gz', fp0, f));
                    end
                end
                sessd = this.sessionData;
                sessd.rnumber = r;
            end            
            delete(sprintf('%s_*_*.4dfp.*', ipr.maskForImages));
        end
        
        %% UTILITY         
              
        function         copySourceToDest(this, ipr)
            if (1 == this.sessionData.rnumber)
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
                this.fileprefixResolved(ipr.dest, this.sessionData.rnumber-1), ...
                this.fileprefixRevision(ipr.dest, this.sessionData.rnumber));
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
        function fp    = fileprefixOfReference(this, ipr)
            fp = this.fileprefixIndexed(mybasename(ipr.dest), this.indexOfReference);
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

 