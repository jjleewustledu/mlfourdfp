classdef CompositeT4ResolveBuilder < mlfourdfp.AbstractT4ResolveBuilder
	%% COMPOSITET4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 21-Jan-2017 13:47:30
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	methods 		          
        function fqfp         = embedInEuclideanR3(this, varargin)
            fqfp = ensureCell( ...                   
                   this.ensureSafeFileprefix(varargin{:}));
            for i = 1:length(fqfp)
                if (isempty(fqfp{i}))
                    continue
                end
                ic = mlfourd.ImagingContext2([fqfp{i} '.4dfp.hdr']);
                if (4 == length(size(ic)) && size(ic,4) > 1) % short-circuit
                    if (lexist([ic.fqfileprefix '_avgt.4dfp.hdr']))
                        ic = mlfourd.ImagingContext2([ic.fqfileprefix '_avgt.4dfp.hdr']);
                    else                        
                        ic = ic.timeAveraged('taus', this.taus);
                        ic.save;
                        fprintf('mlfourdfp.CompositeT4ResolveBuilder.embedInEuclideanR3 saved %s\n', ...
                            ic.fqfilename);
                    end
                end
                fqfp{i} = ic.fqfileprefix;
            end
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
            addParameter(ip, 'logPath',        this.getLogPath,        @ischar);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            this.maskForImages_ = ip.Results.maskForImages;
            this.resolveTag = ip.Results.resolveTag;  
            ipr = ip.Results;        
            ipr = this.expandDest(ipr);
            ipr = this.expandBlurs(ipr);            
            ipr.source = this.ensureLocalFourdfp(ipr.source);
            ipr.resolved = ipr.source; 
            if (~any(this.indicesLogical))
                this = this.copySourceToResolved(ipr);
                this = this.finalize(ipr);
                return
            end
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
            this.imageRegLog = this.loggerFilename( ...
                ipr.dest{this.indexOfReference}, 'func', 'CompositeT4ResolveBuilder_imageReg', 'path', ipr.logPath);
            this.resolveLog = this.loggerFilename( ...
                ipr.dest{this.indexOfReference}, 'func', 'CompositeT4ResolveBuilder_t4ResolveAndPaste', 'path', ipr.logPath);
            
            stagedImgs  = this.lazyStageImages(ipr);    % contracted wrt this.indicesLogical
            if (length(stagedImgs) < 2)                 % degenerate case; proceed to finalize operations
                this.rnumber = this.NRevisions + 1;
                return
            end
            blurredImgs = this.lazyBlurImages(ipr);     % "
            maskedImgs  = this.lazyMasksForImages(ipr); % "
            assertSizeEqual(stagedImgs, blurredImgs, maskedImgs);
            this = this.imageReg(stagedImgs, blurredImgs, maskedImgs);
            [ipr,~,this] = this.resolveAndPaste(ipr); 
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end  
        function this =         imageReg(this, stagedImgs, blurredImgs, maskedImgs)
            len = sum(this.indicesLogical);
            t4fails = zeros(len, len);
            for m = 1:len
                for n = 1:len
                    if (m ~= n)
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            if (~this.valid_t4(t4))
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
                            t4fails(m,n) = t4fails(m,n) + 1;
                            dispwarning(ME, 'mlfourdfp:RuntimeWarning', ...
                                'CompositeT4ResolveBulder.imageReg could not operate %s on dest %s, source %s, t4 %s', ...
                                this.sessionData.compAlignMethod, blurredImgs{n}, blurredImgs{m}, t4);
                        end                        
                    end
                end
            end 
            
            this.t4ResolveError_.logger.add( ...
                sprintf('CompositeT4ResolveBuilder.imageReg.t4fails->\n%s', mat2str(t4fails)));            
            this.indicesLogical(this.indicesLogical) = ...
                ensureRowVector(this.indicesLogical(this.indicesLogical)) & ...
                ensureRowVector(sum(t4fails,1) < 0.25*len);   
            this.t4ResolveError_.logger.add( ...
                sprintf('CompositeT4ResolveBuilder.imageReg.indicesLogical->\n%s', mat2str(this.indicesLogical)));
            [this.t4ResolveError_,this.t4_resolve_err] = ...
                this.t4ResolveError_.estimateErr(stagedImgs, this.indicesLogical, 'rnumber', this.rnumber);
            
            this.deleteTrash;
        end  
        function [ipr,imgFps,this] = resolveAndPaste(this, ipr)
            %% RESOLVEANDPASTE - preassign ipr.dest, this.resolveTag, this.indexOfReference as needed.
            %  @param ipr is a struct w/ field dest, a cell array of fileprefixes || is a cell
            
            if (iscell(ipr))
                ipr_cell = ipr;
                ipr = struct('dest', '');
                ipr.dest = ipr_cell;
            end
            if (~any(this.indicesLogical))
                dest_ = cellfun(@(x) mybasename(x), ipr.dest, 'UniformOutput', false);
                ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest_, 'UniformOutput', false);
                cellfun(@(x,y) this.buildVisitor.copyfile_4dfp(x,y), dest_, ipr.resolved, 'UniformOutput', false);
                imgFps = '';
                this.rnumber = this.NRevisions;
                return
            end            
            
            pwd0   = pushd(fileparts(ipr.dest{1}));
            imgFpsc = {mybasename(this.fileprefixOfReference(ipr))}; % initial on ipr.dest
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    % fileprefix ipr.dest{f} != this.indexOfReference
                    imgFpsc = [imgFpsc mybasename(ipr.dest{f})]; %#ok<AGROW>
                end
            end
            
            %% Must use short fileprefixes in calls to t4_resolve to avoid filenaming error by t4_resolve  
            %  t4_resolve: /data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/V1/FDG_V1-NAC/E8/fdgv1e8r1_frame1_to_/data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/V1/FDG_V1-NAC/E8/fdgv1e8r1_frame8_t4 read error         
            imgFps = cell2str(imgFpsc, 'AsRow', true);
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, imgFps, ...
                'options', '-v -m -s', 'log', this.resolveLog);
            this = this.cacheT4s(imgFpsc); % for each rnumber
            this.t4imgAll(ipr, this.resolveTag); % transform ipr.dest on this.resolveTag
            dest_ = cellfun(@(x) mybasename(x), ipr.dest, 'UniformOutput', false);
            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest_, 'UniformOutput', false); 
            movefile([this.resolveTag '.mat0'], [ipr.resolved{this.indexOfReference} '_' mydatetimestr(now) '.mat0']);
            movefile([this.resolveTag '.sub'],  [ipr.resolved{this.indexOfReference} '_' mydatetimestr(now) '.sub']);
            popd(pwd0);
        end
        function t4           = t4_to_resolveTag(this, varargin)
            %% T4_TO_RESOLVETAG provides the t4 filename to the target of t4_resolve, generating the t4 file as needed.
            %  @param required idx is numeric.
            %  @param named rnumber is numeric.
            %  @param named rstring is char.
            %  @return t4 filename is char.
            
            ip = inputParser;
            addRequired( ip, 'idx', @isnumeric);
            addParameter(ip, 'rnumber', this.NRevisions, @isnumeric);
            addParameter(ip, 'rstring', '', @ischar);
            parse(ip, varargin{:});
            
            t4 = this.getT4__(ip.Results.idx, this.getRstring__(ip.Results));
            if (lexist(t4))
                return
            end
            
            % aufbau composite t4 for all r-numbers
            for r = 1:this.NRevisions
                assert(lexist(this.getT4__(ip.Results.idx, ip.Results.rstring))); 
            end
            rstr = 'r1';
            rstr1 = ip.Results.rstring;
            for r = 1:this.NRevisions-1
                rstr1 = sprintf('%sr%i', rstr, r+1);
                this.buildVisitor.t4_mul( ...
                    this.getT4__(ip.Results.idx, rstr), ...
                    this.getT4__(ip.Results.idx, sprintf('r%i', r+1)), ...
                    this.getT4__(ip.Results.idx, rstr1));
                rstr = rstr1;
            end
            t4 = this.getT4__(ip.Results.idx, rstr1);
        end
        function this         = alreadyFinalized(this, ipr)
            dest         = cellfun(@(x) this.fileprefixRevision(x, this.NRevisions), ipr.dest, 'UniformOutput', false);
            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest, 'UniformOutput', false);
            
            this.logger_.add(sprintf('%s is finished; building product with ipr->\n%s', class(this), evalc('disp(ipr)')));
            this.logger_.add(sprintf('ipr.resolved->\n%s', evalc('disp(ipr.resolved'')')));
            this = this.buildProduct(ipr);
        end
        function                teardownResolve(this, ipr)
            if (this.keepForensics); return; end
            
            for r = 1:this.NRevisions
                for il = 1:length(this.indicesLogical)
                    fp0 = ipr.dest{il};
                    if (this.indicesLogical(il))
                        try
                            deleteExisting(sprintf('%sr%i.4dfp.*',    fp0, r));
                            deleteExisting(sprintf('%sr%i_b*.4dfp.*', fp0, r));
                            deleteExisting(sprintf('%sr%i_C*.4dfp.*', fp0, r));
                            deleteExisting(sprintf('%sr%i_%s.4dfp.*', fp0, r, this.resolveTag));
                            deleteExisting(sprintf('%sr%i_g*.nii.gz', fp0, r));
                        catch ME
                            handwarning(ME, ...
                                'mlfourdfp:RuntimeWarning', ...
                                'CompositeT4ResolveBuilder.teardownResolve failed to delete %sr%i', fp0, r);
                        end
                    end
                end
            end
            if (iscell(ipr.maskForImages))
                cellfun(@(x) deleteExisting(sprintf('%s_*_*.4dfp.*', x)), ipr.maskForImages);
            else
                deleteExisting(sprintf('%s_*_*.4dfp.*', ipr.maskForImages));
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
                                %ipr.dest{s} = sprintf('%s_%s', ipr.dest{s}, mydatetimestr(now));
                            end
                        end
                    end
                catch ME
                    handexcept(ME, 'mlfourdfp:RuntimeError', ...
                        'CompositeT4ResolveBuilder.copySourceToDest could not copy 4dfp %s to %s', ...
                        ipr.source{s}, ipr.dest{s});
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
        function fqfps   = lazyMasksForImages(this, ipr, varargin)
            %  @param ipr.maskForImages is 'none' or the fileprefix/name of an anatomical image which will be 
            %  thresholded to generate a mask; the anatomical image must have transverse orientation.   
            %  See also:  mlfourdfp.FourdfpVisitor.transverse_t4.
            %  @param this.theImages must be populated with the images to be masked.
            %  @return fqfps, a cell-array of binary masks, identically sized to this.theImages.
            
            fqfps = cell(1, length(ipr.source));            
            if (~iscell(ipr.maskForImages) && strcmpi(ipr.maskForImages, 'none'))
                fqfps = cellfun(@(x) 'none', fqfps, 'UniformOutput', false);
                return
            end   
            
            % build masks for each frame
            if (~iscell(ipr.maskForImages))
                ipr.maskForImages = repmat(ensureCell(ipr.maskForImages), size(ipr.source));
            end
            assertSizeEqual(ipr.maskForImages, ipr.source);
            for ii = 1:length(ipr.maskForImages)           
                if (strcmp(ipr.maskForImages{ii}, 'Msktgen'))
                    try
                        mg   = mlpet.Msktgen( ...
                            'sessionData', this.sessionData, ...
                            'logPath', fullfile(pwd, 'Log', ''));
                        mskt = mg.constructMskt( ...
                            'source', ipr.source{ii}, ...
                            'blurArg', ipr.sourceBlur(ii), ...
                            'intermediaryForMask', this.sessionData.T1001, ...
                            'sourceOfMask', fullfile(this.sessionData.sessionPath, 'brainmask.4dfp.hdr'), ...
                            'blurForMask', 33);
                        if (~lexist_4dfp(mskt.fileprefix))
                            this.buildVisitor.copyfile_4dfp(mskt.fqfileprefix, mskt.fileprefix);
                        end
                        fqfps{ii} = mskt.fqfileprefix;
                    catch ME
                        dispexcept(ME, 'mlfourdfp:RuntimeError', ...
                            'CompositeT4ResolveBuilder.lazyMaskForImages failed case Msktgen on %s; trying %s', ...
                            ipr.source{ii}, 'none');
                    end
                    continue
                end
                if (strcmp(ipr.maskForImages{ii}, 'T1001'))
                    try
                        if (~lexist_4dfp([ipr.maskForImages{ii} '_mskt']))
                            this.buildVisitor.msktgenMprage(ipr.maskForImages{ii});
                        end
                        fqfps{ii} = [ipr.maskForImages{ii} '_mskt'];
                    catch ME
                        dispexcept(ME, 'mlfourdfp:RuntimeError', ...
                            'CompositeT4ResolveBuilder.lazyMaskForImages failed case T1001; trying none');
                    end
                    continue
                end
                if (strcmp(ipr.maskForImages{ii}, 'msktgen_4dfp'))
                    try
                        if (~lexist_4dfp([ipr.source{ii} '_mskt']))
                            this.buildVisitor.msktgenMprage(ipr.source{ii}, this.sessionData.atlas('typ','fqfp'));
                        end
                        fqfps{ii} = [ipr.source{ii} '_mskt'];
                    catch ME
                        dispexcept(ME, 'mlfourdfp:RuntimeError', ...
                            'CompositeT4ResolveBuilder.lazyMaskForImages failed case msktgen_4dfp on %s; trying %s', ...
                            ipr.source{ii}, 'none');
                    end
                    continue
                end
                if (~isempty(ipr.maskForImages{ii}))
                    fqfps{ii} = ipr.maskForImages{ii};
                    continue
                end
                fqfps{ii} = 'none';
            end
        end  
        function        cd(~, varargin)
            if (isempty(varargin) || isempty(varargin{1}))
                return
            end
            cd(varargin{:});
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
        
        %%
        
 		function this = CompositeT4ResolveBuilder(varargin)
 			%% COMPOSITET4RESOLVEBUILDER
            %  TODO:  consider support for cctor;  See also mlfourdfp.T4ResolveBuilder.
            %  @param theImages =: safe fileprefix =: time summed.  Always contracted to ${\rm I\!R}^3$ by method
            %  embedInEuclideanR3.  
            %  @param blurArg; default := this.sessiondata.compositeT4ResolveBuilderBlurArg.
            %  @param indicesLogical 
            %  @param indexOfReference; default := 1, the first of theImages.
 			
            this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'blurArg', this.sessionData.compositeT4ResolveBuilderBlurArg, @isnumeric);
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});
                                        
            this.imageComposite_ = mlfourdfp.ImageComposite( ...
                this, ...
                'theImages', this.embedInEuclideanR3( ...
                             this.ensureLocalFourdfp(this.theImages)), ...
                'indicesLogical', ip.Results.indicesLogical, ...
                'indexOfReference', ip.Results.indexOfReference);
            this.blurArg_ = ip.Results.blurArg;  
            this = this.updateFinished;
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
                    this.product_{il1} = mlfourd.ImagingContext2([ipr.resolved{il1} '.4dfp.hdr']);
                    assert(~isempty(ipr.resolved{il1}), ...
                        'mlfourdfp:ValueError', ...
                        'T4ResolveBuilder.buildProduct.ipr.resolved{%i}->%s', il1, ipr.resolved{il1});
                    assert(~isempty(this.product_{il1}), ...
                        'mlfourdfp:ValueError', ...
                        'T4ResolveBuilder.buildProduct.this.product{%i}->%s', il1, this.product_{il1}.fqfilename);
                    this.product_{il1}.addImgrec( ...
                        ['mlfourdfp.CompositeT4ResolveBuilder.buildProduct()\nipr->\n' evalc('disp(ipr)')]);
                    this.product_{il1}.addLog( ...
                        ['mlfourdfp.CompositeT4ResolveBuilder.t4ResolveError_.logger.fqfilename->' ...
                        this.t4ResolveError_.logger.fqfilename]);
                end
            end            
        end
        function this = cacheT4s(this, imgFpsc)
            %  @return this.t4s_{r}{1} for r-number r is the reference; size(this.t4s_) == this.NRevisions;
            %  size(this.t4s_{r}) == size(this.indicesLogical).
            
            this.t4s_{this.rnumber} = cell(size(imgFpsc));
            for f = 1:length(imgFpsc)
                try
                    this.t4s_{this.rnumber}{f} = sprintf('%s_to_%s_t4', imgFpsc{f}, this.resolveTag);
                catch ME
                    dispwarning(ME, 'mlfourdfp:RuntimeWarning', ...
                        'AbstractT4ResolveBuilder.cacheT4s could not assign %s', this.t4s_{this.rnumber}{f});
                end                
            end
        end
    end

    %% PRIVATE
    
    methods (Access = private)
        function rstr_ = getRstring__(this, ipr_)
            if (isempty(ipr_.rstring))
                rstr_ = 'r1';
                for r_ = 1:this.NRevisions-1
                    rstr_ = sprintf('%sr%i', rstr_, r_+1);
                end
            else
                rstr_ = ipr_.rstring;
            end
        end
        function t4_   = getT4__(this, idx_, rstr_)
            t4_ = sprintf('%s%s_to_%s_t4', ...
                this.imageComposite.theImages{idx_}, rstr_, this.resolveTag);
        end
    end
        
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

