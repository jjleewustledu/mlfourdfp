classdef CompositeT4ResolveBuilder < mlfourdfp.AbstractT4ResolveBuilder
	%% COMPOSITET4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 21-Jan-2017 13:47:30
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	methods 		  
 		function this = CompositeT4ResolveBuilder(varargin)
 			%% COMPOSITET4RESOLVEBUILDER
            %  @param blurArg; default := this.sessiondata.compositeT4ResolveBuilderBlurArg.
            %  @param theImages =: safe fileprefix =: time summed.
            %  @param indicesLogical 
            %  @param indexOfReference; default := 1, the first of theImages.
            %  @param maskForImages is a fileprefix or 'none'.
 			
            this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'blurArg', this.sessionData.compositeT4ResolveBuilderBlurArg, @isnumeric);
            addParameter(ip, 'theImages', {}, @(x) ~isempty(x));
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            addParameter(ip, 'maskForImages', 'none', @(x) ~isempty(x));
            parse(ip, varargin{:});
            
            this.imageComposite_ = mlfourdfp.ImageComposite( ...
                this, ...
                'theImages', this.ensureSafeFileprefix( ...
                             this.embedInEuclideanR3(ip.Results.theImages)), ...
                'indicesLogical', ip.Results.indicesLogical, ...
                'indexOfReference', ip.Results.indexOfReference);
            this.indexOfReference = ip.Results.indexOfReference;
            this.blurArg_ = ip.Results.blurArg;  
            this.maskForImages_ = ip.Results.maskForImages;
            this = this.updateFinished;
        end
        function imgs         = embedInEuclideanR3(this, varargin)
            imgs = ensureCell(varargin);
        end
        
        function this         = resolve(this, varargin)
            import mlfourdfp.*;
            ip = inputParser;
            addParameter(ip, 'dest',           {},                  @iscell);
            addParameter(ip, 'source',         this.theImages,      @iscell);
            addParameter(ip, 'destMask',       'none',              @ischar);
            addParameter(ip, 'sourceMask',     'none',              @ischar);
            addParameter(ip, 'destBlur',       this.blurArg,        @isnumeric); % fwhh/mm
            addParameter(ip, 'sourceBlur',     this.blurArg,        @isnumeric); % fwhh/mm
            addParameter(ip, 'maskForImages',  this.maskForImages_, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 't40',            this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'resolveTag',     this.resolveTag,     @ischar);
            addParameter(ip, 'log',            '/dev/null',         @ischar);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            this.resolveTag = ip.Results.resolveTag;  
            ipr = ip.Results;            
            ipr = this.expandDest(ipr);
            ipr = this.expandBlurs(ipr);            
            ipr.resolved = ipr.source; 
            if (this.isfinished)  
                this = this.alreadyFinalized(ipr);
                return
            end
            while (this.rnumber <= this.NRevisions)
                ipr.source = ipr.resolved;
                ipr.dest   = cellfun(@(x) this.fileprefixRevision(x, this.rnumber), ipr.dest, 'UniformOutput', false);
                [ipr,this] = this.revise(ipr);
                assert(this.rnumber < 10);
            end  
            this = this.finalize(ipr);
        end        
        function [ipr,this]   = revise(this, ipr)
            ipr = this.copySourceToDest(ipr); % crop/copy ipr.source to ipr.dest       
            this.imageRegLog = loggerFilename( ...
                ipr.dest{this.indexOfReference}, 'func', 'CompositeT4ResolveBuilder_imageReg', 'path', this.logPath);
            this.resolveLog = loggerFilename( ...
                ipr.dest{this.indexOfReference}, 'func', 'CompositeT4ResolveBuilder_t4ResolveAndPaste', 'path', this.logPath);
            
            this = this.imageReg(ipr);
            ipr = this.resolveAndPaste(ipr); 
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end  
        function this =         imageReg(this, ipr)
            stagedImgs   = this.lazyStageImages(ipr);
            blurredImgs  = this.lazyBlurImages(ipr);
            maskedImgs     = this.lazyMasksForImages(ipr);
            assertSizeEqual(stagedImgs, blurredImgs, maskedImgs);
            len = length(stagedImgs);
            t4Failures = zeros(len, len);
            for m = 1:len
                for n = 1:len
                    if (m ~= n && ...
                        this.indicesLogical(m) && this.indicesLogical(n)) 
                    
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            if (~lexist(t4))
                                this.buildVisitor.(this.sessionData.compAlignMethod)( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
                                    'destMask',   maskedImgs{m}, ...
                                    'sourceMask', maskedImgs{n}, ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'log',        this.imageRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of image files
                            % e. g., image1_to_image2_t4    
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
            fprintf('CompositeT4ResolveBuilder.imageReg: t4Failures->%s\n', mat2str(t4Failures));
            this.indicesLogical(this.indicesLogical) = ...
                ensureRowVector(this.indicesLogical(this.indicesLogical)) & ...
                ensureRowVector(t4Failures < 0.25*len);   
            fprintf('CompositeT4ResolveBuilder.imageReg: this.indicesLogical->%s\n', mat2str(this.indicesLogical)); 
            
            this.t4_resolve_err = nan(len, len);
            for m = 1:length(stagedImgs)
                for n = 1:length(stagedImgs)
                    if (m ~= n && ...
                        this.indicesLogical(m) && this.indicesLogical(n)) 
                    
                        try               
                            [rmsdeg,rmsmm] = this.t4_resolve_errParser(this.resolvePair( ...
                                    mybasename(ipr.dest{m}), mybasename(ipr.dest{n})));
                            this.t4_resolve_err(m,n) = this.t4_resolve_errAverage(rmsdeg, rmsmm);
                        catch ME
                            dispwarning(ME);
                        end
                        
                    end
                end
            end 
            fprintf('CompositeT4ResolveBuilder.imageReg: this.t4_resolve_err->%s\n', mat2str(this.t4_resolve_err));
            
            this.deleteTrash;
        end  
        function r            = resolvePair(this, f1, f2)
            [~,r] = this.buildVisitor.t4_resolve( ...
                this.resolveTag, [f1 ' ' f2], ...
                'options', '-v -m -s');
        end
        function [ipr,imgFps] = resolveAndPaste(this, ipr)
            %% RESOLVEANDPASTE - preassign ipr.dest, this.resolveTag, this.indexOfReference as needed.
            %  @param ipr is a struct w/ field dest, a cell array of fileprefixes || is a cell
            
            if (iscell(ipr))
                ipr_cell = ipr;
                ipr = struct('dest', '');
                ipr.dest = ipr_cell;
            end
            pwd0   = pushd(fileparts(ipr.dest{1}));
            imgFps = mybasename(this.fileprefixOfReference(ipr)); % initial on ipr.dest
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    %                    fileprefix of frame != this.indexOfReference
                    imgFps = [imgFps ' ' mybasename(ipr.dest{f})]; %#ok<AGROW>
                end
            end
                     
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, imgFps, ...
                'options', '-v -m -s', 'log', this.resolveLog);
            this.t4imgAll(ipr, this.resolveTag); % transform ipr.dest on this.resolveTag
            dest_ = cellfun(@(x) mybasename(x), ipr.dest, 'UniformOutput', false);
            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest_, 'UniformOutput', false); 
            movefile([this.resolveTag '.mat0'], [ipr.resolved{this.indexOfReference} '_' datestr(now, 30) '.mat0']);
            movefile([this.resolveTag '.sub'],  [ipr.resolved{this.indexOfReference} '_' datestr(now, 30) '.sub']);
            popd(pwd0);
        end
        function t4           = t4_to_resolveTag(this, varargin)
            ip = inputParser;
            addRequired( ip, 'idx', @isnumeric);
            addParameter(ip, 'rnumber', this.NRevisions, @isnumeric);
            addParameter(ip, 'rstring', '', @ischar);
            parse(ip, varargin{:});
            
            t4 = getT4__(ip.Results.idx, getRstring__(ip.Results));
            if (lexist(t4))
                return
            end
            
            % aufbau composite t4 for all r-numbers
            for r = 1:this.NRevisions
                assert(lexist(this.t4_to_resolveTag(ip.Results.idx, 'rnumber', r)));
            end
            rstr = 'r1';
            for r = 1:this.NRevisions-1
                rstr1 = sprintf('%sr%i', rstr, r+1);
                this.buildVisitor.t4_mul( ...
                    getT4__(this, ip.Results.idx, rstr), ...
                    getT4__(this, ip.Restuls.idx, sprintf('r%i', r+1)), ...
                    getT4__(this, ip.Results.idx, rstr1));
                rstr = rstr1;
            end
            t4 = getT4__(this, ip.Results.idx, rstr1);
            
            function rstr_ = getRstring__(ipr_)
                if (isempty(ipr_.rstring))
                    rstr_ = sprintf('r%i', ipr_.rnumber);
                else
                    rstr_ = ipr_.rstring;
                end
            end
            function t4_   = getT4__(this, idx_, rstr_)
                t4_ = sprintf('%s%s_to_%s_t4', ...
                    this.imageComposite.theImages{idx_}, rstr_, this.resolveTag);
            end
        end
        function this         = alreadyFinalized(this, ipr)
            dest         = cellfun(@(x) this.fileprefixRevision(x, this.NRevisions), ipr.dest, 'UniformOutput', false);
            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest, 'UniformOutput', false);
            this         = this.buildProduct(ipr);
        end
        function                teardownResolve(this, ipr)
            if (this.keepForensics); return; end
            
            try
                for r = 1:this.NRevisions                
                    for il = 1:length(this.indicesLogical)
                        fp0 = ipr.dest{il};
                        if (this.indicesLogical(il))
                            deleteExisting(sprintf('%sr%i.4dfp.*', fp0, r));
                            deleteExisting(sprintf('%sr%i_b*.4dfp.*', fp0, r));
                            deleteExisting(sprintf('%sr%i_C*.4dfp.*', fp0, r));
                            deleteExisting(sprintf('%sr%i_%s.4dfp.*', fp0, r, this.resolveTag));
                            deleteExisting(sprintf('%sr%i_g*.nii.gz', fp0, r));
                        end
                    end
                end            
                deleteExisting(sprintf('%s_*_*.4dfp.*', ipr.maskForImages));
            catch ME
                handwarning(ME)
            end
        end
        
        %% UTILITY
        
        function ipr     = copySourceToDest(this, ipr)
            if (this.skipT4imgAll)
                return
            end
            
            if (1 == this.rnumber)
                try              
                    for s = 1:length(ipr.source)
                        if (this.indicesLogical(s))
                            ipr.dest{s} = this.fileprefixRevision(ipr.source{s}, 1);
                            if (~this.buildVisitor.lexist_4dfp(ipr.dest{s}) && ...
                                ~strcmp(ipr.dest{s}, ipr.source{s}))
                                this.buildVisitor.copy_4dfp(ipr.source{s}, ipr.dest{s});
                                %ipr.dest{s} = sprintf('%s_%s', ipr.dest{s}, datestr(now, 30));
                            end
                        end
                    end
                catch ME
                    handexcept(ME);
                end
                return
            end
            this.copyResolvedToNewRevision(ipr);
        end
        function           copyResolvedToNewRevision(this, ipr)
            %% COPYRESOLVEDTONEWREVISION opportunistically reuses existing files from the last iteration
            
            for s = 1:length(ipr.source)
                if (this.indicesLogical(s))
                    this.buildVisitor.copy_4dfp( ...
                        this.fileprefixResolved(ipr.dest{s}, this.rnumber-1), ...
                        this.fileprefixRevision(ipr.dest{s}, this.rnumber));
                end
            end
        end
        function ipr     = expandDest(this, ipr)
            if (isempty(ipr.dest))
                for s = 1:length(ipr.source)
                    ipr.dest{s} = sprintf('%s_%s', ipr.source{s}, this.resolveTag);
                end
            end
            assert(length(ipr.source) == length(ipr.source));
        end
        function fp      = fileprefixIndexed(~, ipr)
            assert(isfield(ipr, 'dest'));
            assert(isfield(ipr, 'currentIndex'));
            fp = ipr.dest{ipr.currentIndex};
        end
        function fp      = fileprefixOfReference(this, ipr)
            fp = ipr.dest{this.indexOfReference};
        end
        function lazyMsk = lazyMasksForImages(this, ipr)
            %  @param ipr.maskForImages is 'none' or the fileprefix/name of an anatomical image which will be 
            %  thresholded to generate a mask; the anatomical image must have transverse orientation.   
            %  See also:  mlfourdfp.FourdfpVisitor.transverse_t4.
            %  @param this.theImages must be populated with the images to be masked.
            %  @return mska, a cell-array of binary masks, identically sized to this.theImages.
            
            lazyMsk = cell(1, length(ipr.source));
            
            if (strcmp(ipr.maskForImages, 'none'))
                lazyMsk = cellfun(@(x) 'none', lazyMsk, 'UniformOutput', false);
                return
            end
            
            assert(iscell(ipr.maskForImages));
            assert(length(ipr.maskForImages) == length(ipr.source));
            for ii = 1:length(ipr.source)
                if (strcmp(ipr.maskForImages{ii}, 'none'))
                    lazyMsk{ii} = 'none';
                    continue
                end
                if (strcmp(ipr.maskForImages{ii}, 'T1001'))
                    if (~lexist_4dfp( [ipr.maskForImages{ii} '_mskt']))
                        this.buildVisitor.msktgenMprage(ipr.maskForImages{ii});
                    end
                    lazyMsk{ii} = [ipr.maskForImages{ii} '_mskt'];
                    continue
                end
                lazyMsk{ii} = 'none';
            end
        end      
        function fqfps   = lazyStageImages(this, ipr)
            assert(iscell(ipr.dest));
            fqfps = {};
            for d = 1:length(ipr.dest)
                if (this.indicesLogical(d))
                    fqfps = [fqfps ipr.dest{d}]; %#ok<AGROW>
                end
            end
        end        
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function this = buildProduct(this, ipr)
            this.ipResults_ = ipr;
            this.rnumber = this.NRevisions;            
            this.product_ = cell(1, sum(this.indicesLogical));            
            il1 = 0;
            for il = 1:length(this.indicesLogical)
                if (this.indicesLogical(il))
                    il1 = il1 + 1;                    
                    this.product_{il1} = mlfourd.ImagingContext([ipr.resolved{il} '.4dfp.ifh']);
                end
            end            
        end
    end

    %% PRIVATE

    properties (Access = private)
        maskForImages_
    end
        
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

