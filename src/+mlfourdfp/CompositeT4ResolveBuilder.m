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
 			
            this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'blurArg', this.sessionData.compositeT4ResolveBuilderBlurArg, @isnumeric);
            addParameter(ip, 'theImages', {}, @(x) ~isempty(x));
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            addParameter(ip, 'maskForImages', 'none', @ischar);
            parse(ip, varargin{:});
            
            import mlfourdfp.*;
            this.imageComposite_ = ImageComposite(this, ...
                'theImages', FourdfpVisitor.ensureSafeFileprefix(ensureCell(ip.Results.theImages)), ...
                'indicesLogical', ip.Results.indicesLogical, ...
                'indexOfReference', ip.Results.indexOfReference);
            this.indexOfReference = ip.Results.indexOfReference;
            this.blurArg_ = ip.Results.blurArg;  
            this.maskForImages_ = ip.Results.maskForImages;
            this = this.updateFinished;
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
            addParameter(ip, 'maskForImages',  this.maskForImages_, @ischar);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 't40',            this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'resolveTag',     this.resolveTag,     @ischar);
            addParameter(ip, 'log',            '/dev/null',         @ischar);
            parse(ip, varargin{:});
            %if (this.isfinished); return; end
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
            
            this.imageReg(ipr);
            ipr = this.resolveAndPaste(ipr); 
            this.teardownRevision(ipr);
            this.rnumber = this.rnumber + 1;
        end  
        function                imageReg(this, ipr)
            stagedImgs   = this.lazyStageImages(ipr);
            blurredImgs  = this.lazyBlurImages(ipr);
            maskImgs     = this.lazyMasksForImages(ipr);
            assert(length(stagedImgs) == length(blurredImgs));
            for m = 1:length(stagedImgs)
                for n = 1:length(stagedImgs)
                    if (m ~= n) 
                        try
                            t4 = this.buildVisitor.filenameT4(stagedImgs{n}, stagedImgs{m});
                            if (~lexist(t4))
                                this.buildVisitor.(this.sessionData.compAlignMethod)( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
                                    'destMask',   maskImgs{m}, ...
                                    'sourceMask', maskImgs{n}, ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'log',        this.imageRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of image files
                            % e. g., image1_to_image2_t4                            
                        catch ME
                            handerror(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
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
        function dest1        = reconstituteImages(this, ipr, varargin)
            if (this.skipT4imgAll)
                return
            end
            
            ip = inputParser;
            addRequired(ip, 'ipr', @(x) isstruct(x) && ...
                                        isfield(x, 'fqfps') && length(x.fqfps) > 1 && ...
                                        isfield(x, 'dest1') && ...
                                        ischar(x.dest));
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, ipr, varargin{:});

            pasteList = [ipr.dest1 '_paste.lst'];
            if (lexist(pasteList))
                delete(pasteList); 
            end
            fid = fopen(pasteList, 'w');
            for idx = 1:length(ipr.fqfps)
                fprintf(fid, '%s.4dfp.img\n', ipr.fqfps{idx});
            end
            fclose(fid);            
            this.buildVisitor.paste_4dfp(pasteList, ipr.dest1, 'options', '-a ');
            dest1 = ipr.dest1;
            %deletefiles(fqfps{:});
        end  
        function this         = finalize(this, ipr)
            this = this.buildProduct(ipr);
            this.teardownResolve(ipr);
            this.finished.touchFinishedMarker;  
        end
        function this         = alreadyFinalized(this, ipr)    
            dest = cellfun(@(x) this.fileprefixRevision(x, this.NRevisions), ipr.dest, 'UniformOutput', false);
            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest, 'UniformOutput', false);
            this = this.buildProduct(ipr);
        end
        function this         = buildProduct(this, ipr)      
            this.ipResults_ = ipr;
            this.rnumber = this.NRevisions;            
            this.product_ = cell(1, sum(this.indicesLogical));            
            il1 = 0;
            for il = 1:length(this.indicesLogical)
                if (this.indicesLogical(il))
                    il1 = il1 + 1;                    
                    this.product_{il1} = mlpet.PETImagingContext([ipr.resolved{il} '.4dfp.ifh']);
                    %assert(~isempty(this.product_{il1}));
                end
            end            
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
        
        function ipr   = copySourceToDest(this, ipr)
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
        function         copyResolvedToNewRevision(this, ipr)
            %% COPYRESOLVEDTONEWREVISION opportunistically reuses existing files from the last iteration
            
            for s = 1:length(ipr.source)
                if (this.indicesLogical(s))
                    this.buildVisitor.copy_4dfp( ...
                        this.fileprefixResolved(ipr.dest{s}, this.rnumber-1), ...
                        this.fileprefixRevision(ipr.dest{s}, this.rnumber));
                end
            end
        end
        function ipr   = expandDest(this, ipr)
            if (isempty(ipr.dest))
                for s = 1:length(ipr.source)
                    ipr.dest{s} = sprintf('%s_%s', ipr.source{s}, this.resolveTag);
                end
            end
            assert(length(ipr.source) == length(ipr.source));
        end
        function fp    = fileprefixOfReference(this, ipr)
            fp = ipr.dest{this.indexOfReference};
        end
        function fqfp  = lazyBlurImage(this, ipr, blur)
            %% LAZYBLURIMAGE uses specifiers in ipr; will not replace any existing image
            %  @param ipr is a struct
            %  @param blur is numeric
            
            fqfp_ = ipr.dest{ipr.currentIndex};
            fqfp  = this.fileprefixBlurred(fqfp_, blur);
            if (~this.buildVisitor.lexist_4dfp(fqfp))
                this.buildVisitor.imgblur_4dfp(fqfp_, blur);
            end
        end
        function mska  = lazyMasksForImages(this, ipr)
            %  @param ipr.maskForImages is 'none' or the fileprefix/name of an anatomical image which will be 
            %  thresholded to generate a mask; the anatomical image must have transverse orientation.   
            %  See also:  mlfourdfp.FourdfpVisitor.transverse_t4.
            %  @param this.theImages must be populated with the images to be masked.
            %  @return mska, a cell-array of binary masks, identically sized to this.theImages.
            
            if (strcmp(ipr.maskForImages, 'none'))
                mska = cell(1, length(ipr.dest));
                for f = 1:length(mska)
                    mska{f} = 'none';
                end
                return
            end
            
            msk0    = myfileprefix(ipr.maskForImages);
            [~,msk] = fileparts(this.maskFromImage(msk0));
            bv      = this.buildVisitor;
            theDest = bv.ensureSafeFileprefix(ensureCell(ipr.dest));
            mska = cell(1, length(theDest));
            for ii = 1:length(theDest)
                mska{ii} = ['mask_' theDest{ii}];
                if (~bv.lexist_4dfp(mska{ii}))
                    t4 = bv.align_multiSpectral( ...
                        'dest',       theDest{ii}, ...
                        'source',     msk0, ...
                        'destBlur',   this.blurArg, ...
                        'sourceBlur', this.blurArg, ...
                        't4img_4dfp', false, ...
                        'log',        this.imageRegLog);
                    bv.t4img_4dfp(t4, msk, 'out', mska{ii}, 'options', ['-O' theDest{ii}]);
                end
            end
        end
        function fqfps = lazyStageImages(this, ipr)
            assert(iscell(ipr.dest));
            fqfps = {};
            for d = 1:length(ipr.dest)
                if (this.indicesLogical(d))
                    fqfps = [fqfps ipr.dest{d}]; %#ok<AGROW>
                end
            end
        end
        
        function this  = prepareMprToAtlasT4(this)
            %% PREPAREMPRTOATLAST4
            %  @param this.sessionData.{mprage,atlas} are valid.
            %  @return this.product_ := [mprage '_to_' atlas '_t4'], existing in the same folder as mprage.
            
            sessd      = this.sessionData;
            mpr        = sessd.mprage('typ', 'fp');
            mprToAtlT4 = [mpr '_to_' sessd.atlas('typ', 'fp') '_t4'];            
            if (~lexist(fullfile(sessd.mprage('typ', 'path'), mprToAtlT4)))
                pwd0 = pushd(sessd.mprage('typ', 'path'));
                this.msktgenMprage(mpr);
                popd(pwd0);
            end
            this.product_ = mprToAtlT4;
        end        
    end 

    %% PRIVATE

    properties (Access = private)
        maskForImages_
    end
        
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

