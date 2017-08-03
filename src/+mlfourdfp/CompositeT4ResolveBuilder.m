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
            addParameter(ip, 'blurArg', 1.5, @isnumeric);
            addParameter(ip, 'theImages', {}, @(x) ~isempty(x));
            addParameter(ip, 'indicesLogical', true, @islogical);
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});
            
            import mlfourdfp.*;
            this.imageComposite_ = ImageComposite(this, ...
                'theImages', FourdfpVisitor.ensureSafeFileprefix(ensureCell(ip.Results.theImages)), ...
                'indicesLogical', ip.Results.indicesLogical, ...
                'indexOfReference', ip.Results.indexOfReference);
            this.indexOfReference = ip.Results.indexOfReference;
            this.blurArg_ = ip.Results.blurArg; 
            this.finished = mlpipeline.Finished(this, 'path', this.logPath, 'tag', lower(this.sessionData.tracer));
            %cd(this.sessionData.tracerLocation);
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
            addParameter(ip, 'dest',           {},                  @iscell);
            addParameter(ip, 'source',         this.theImages,      @iscell);
            addParameter(ip, 'destMask',       'none',              @ischar);
            addParameter(ip, 'sourceMask',     'none',              @ischar);
            addParameter(ip, 'destBlur',       this.blurArg,        @isnumeric); % fwhh/mm
            addParameter(ip, 'sourceBlur',     this.blurArg,        @isnumeric); % fwhh/mm
            addParameter(ip, 'maskForImages',  'maskForImages',     @ischar);
            addParameter(ip, 'indicesLogical', this.indicesLogical, @islogical);
            addParameter(ip, 't40',            this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            %addParameter(ip, 'resolveTag',     this.resolveTag,     @ischar);
            addParameter(ip, 'log',            '/dev/null',         @ischar);
            parse(ip, varargin{:});
            if (this.isfinished); return; end
            this.indicesLogical = ip.Results.indicesLogical;
            %this.resolveTag = ip.Results.resolveTag;  
            ipr = ip.Results;            
            ipr = this.expandDest(ipr);
            ipr = this.expandBlurs(ipr);
            
            ipr.resolved = ipr.source;
            while (this.sessionData.rnumber <= this.NRevisions)
                ipr.source = ipr.resolved;
                ipr.dest   = cellfun(@(x) this.fileprefixRevision(x, this.sessionData.rnumber), ipr.dest, 'UniformOutput', false);
                [ipr,this] = this.revise(ipr);
                assert(this.sessionData.rnumber < 10);
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
                                this.buildVisitor.align_multiSpectral( ...
                                    'dest',       blurredImgs{m}, ...
                                    'source',     blurredImgs{n}, ...
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
        function [ipr,imgFns] = resolveAndPaste(this, ipr)
            %% RESOLVEANDPASTE
            %  @param ipr is a struct w/ field dest, a cell array of fileprefixes || is a cell
            
            if (iscell(ipr))
                ipr_cell = ipr;
                ipr = struct('dest', '');
                ipr.dest = ipr_cell;
            end
            dest_ = cellfun(@(x) mybasename(x), ipr.dest, 'UniformOutput', false);
            imgFns = this.fileprefixOfReference(ipr); % initial 
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f) && f ~= this.indexOfReference)
                    imgFns = [imgFns ' ' dest_{f}]; %#ok<AGROW>
                end
            end
                     
            this.buildVisitor.t4_resolve( ...
                this.resolveTag, imgFns, ...
                'options', '-v -m -s', 'log', this.resolveLog);
            this.t4imgAll(ipr, this.resolveTag);

            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest_, 'UniformOutput', false); 
            movefile([this.resolveTag '.mat0'], [ipr.resolved{this.indexOfReference} '_' datestr(now, 30) '.mat0']);
            movefile([this.resolveTag '.sub'],  [ipr.resolved{this.indexOfReference} '_' datestr(now, 30) '.sub']);
        end
        function dest1        = reconstituteImages(this, ipr, varargin)
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
            this.product_ = cell(1, sum(this.indicesLogical));
            il1 = 0;
            for il = 1:length(this.indicesLogical)
                if (this.indicesLogical(il))
                    il1 = il1 + 1;
                    this.product_{il1} = mlpet.PETImagingContext([ipr.resolved{il} '.4dfp.ifh']);
                    this.buildVisitor.imgblur_4dfp(ipr.resolved{il}, this.blurArg)
                end
            end
            this.teardownResolve(ipr);
            this.finished.touchFinishedMarker;  
        end
        function                teardownResolve(this, ipr)
            if (this.keepForensics); return; end
            
            for r = 1:this.NRevisions                
                for il = 1:length(this.indicesLogical)
                    fp0 = ipr.dest{il};
                    if (this.indicesLogical(il))
                        %delete(sprintf('%sr%i.4dfp.*', fp0, r));
                        delete(sprintf('%sr%i_b*.4dfp.*', fp0, r));
                        delete(sprintf('%sr%i_C*.4dfp.*', fp0, r));
                        %delete(sprintf('%sr%i_%s.4dfp.*', fp0, r, this.resolveTag));
                        delete(sprintf('%sr%i_g*.nii.gz', fp0, r));
                    end
                end
                sessd = this.sessionData;
                sessd.rnumber = r;
            end            
            delete(sprintf('%s_*_*.4dfp.*', ipr.maskForImages));
        end
        
        %% UTILITY
        
        function ipr   = copySourceToDest(this, ipr)
            if (1 == this.sessionData.rnumber)
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
                        this.fileprefixResolved(ipr.dest{s}, this.sessionData.rnumber-1), ...
                        this.fileprefixRevision(ipr.dest{s}, this.sessionData.rnumber));
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
        function fp    = lazyMaskForImages(this, maskFp, fpm, fpn, ~, ~)
            if (isempty(maskFp))
                fp = '';
                return
            end
            
            maskMN = sprintf('%s_%s_%s.4dfp.img', maskFp, fpm, fpn);
            maskNM = sprintf('%s_%s_%s.4dfp.img', maskFp, fpn, fpm);
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
        function fqt4  = mapcycle( this, indexMin, indexMax, varargin)
            switch (this.NRevisions)
                case 1
                    fqt4 = this.map1cycle(indexMin, indexMax, varargin{:});
                case 2
                    fqt4 = this.map2cycle(indexMin, indexMax, varargin{:});
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', ...
                          'HeterogeneousT4ResolveBuilder.mapcycle.NRevisions -> %i', ...
                          this.NRevisions);
            end
            assert(lexist(fqt4, 'file'));
        end
        function fqt4  = map1cycle(this, indexMin, indexMax, varargin)
            ip = inputParser;
            addRequired(ip, 'indexMin', @ischar);
            addRequired(ip, 'indexMax', @ischar);
            addOptional(ip, 't4output', sprintf('map1cycle_to_%sr1_t4', indexMax), @ischar);
            parse(ip, indexMin, indexMax, varargin{:});
            
            pwd0 = pushd(this.sessionData.tracerT4Location);
            b    = this.buildVisitor;
            fqt4 = this.buildVisitor.t4_mul( ...
                     this.T(indexMin, 1), ...
                     b.t4_inv(this.T(indexMax, 1)), basename(ip.Results.t4output));
            fqt4 = fullfile(this.sessionData.tracerT4Location, fqt4);
            popd(pwd0);
        end
        function fqt4  = map2cycle(this, indexMin, indexMax, varargin)
            ip = inputParser;
            addRequired(ip, 'indexMin', @ischar);
            addRequired(ip, 'indexMax', @ischar);
            addOptional(ip, 't4output', sprintf('map2cycle_to_%sr2_t4',indexMax), @ischar);
            parse(ip, indexMin, indexMax, varargin{:});
            
            b = this.buildVisitor;
            
            % 0 := indexMin, F := indexMax, r := resolved
            pwd0 = pushd(this.sessionData.tracerT4Location);
            T_0r_1 = this.T(indexMin, 1);
            T_0r_2 = this.T(indexMin, 2);
            T_rF_2 = b.t4_inv(this.T(indexMax, 2));
            T_rF_1 = b.t4_inv(this.T(indexMax, 1));
            fqt4   = b.t4_mul(T_0r_1, b.t4_mul(T_0r_2, b.t4_mul(T_rF_2, T_rF_1)), basename(ip.Results.t4output));  
            fqt4   = fullfile(this.sessionData.tracerT4Location, fqt4);
            popd(pwd0);
        end
        function t4    = T(this, varargin)
            %% T selects an affine transformation file.
            %  @param f1 is the originating frame/image (char).
            %  @param r1 is the originating iteration of t4_resolve (numeric).
            %  @param f2 is the destination frame/image (char).
            %  @param r2 is the destination iteration of t4_resolve (numeric).
            %  @param fq specifies whether to return the transformation file fully-qualified (logical).
            
            ip = inputParser;
            addRequired( ip, 'f1', @ischar);
            addRequired( ip, 'r1', @isnumeric);
            addOptional( ip, 'f2', '',    @ischar);
            addOptional( ip, 'r2', nan,   @isnumeric);
            addParameter(ip, 'fq', false, @islogical);
            parse(ip, varargin{:});
            
            if (isempty(ip.Results.f2) || isnan(ip.Results.r2))
                t4__ = sprintf('%sr%i_to_resolved_t4', ...
                    ip.Results.f1, ip.Results.r1);
            else
                t4__ = sprintf('%sr%i_to_%sr%i_t4', ...
                    ip.Results.f1, ip.Results.r1, ip.Results.f2, ip.Results.r2);
            end    
            if (~ip.Results.fq)
                t4 = t4__;
            else
                t4 = fullfile(this.sessionData.tracerT4Location, t4__);                
            end
        end                  
        function this  = t4imgAll(this, ipr, tag)
            tag = mybasename(tag);
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    destFp = ipr.dest{f};
                    this.buildVisitor.t4img_4dfp( ...
                        sprintf('%s_to_%s_t4', destFp, tag), ...
                        destFp, ...
                        'out', sprintf('%s_%s', destFp, tag), ...
                        'options', ['-O' ipr.dest{this.indexOfReference}]);
                end
            end
        end  
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

