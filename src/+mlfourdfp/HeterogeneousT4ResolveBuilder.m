classdef HeterogeneousT4ResolveBuilder < mlfourdfp.T4ResolveBuilder
    %% HETEROGENEOUST4RESOLVEBUILDER
    %  @deprecated prefer CompositeT4ResolveBuilder.
    %  JJL.  2017 Aug 4.
    
    
	methods 
		  
 		function this = HeterogeneousT4ResolveBuilder(varargin)
 			%% HETEROGENEOUST4RESOLVEBUILDER
 			
            this = this@mlfourdfp.T4ResolveBuilder(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'blurArg', 2.03, @isnumeric);
            addParameter(ip, 'indicesLogical', [], @islogical);
            parse(ip, varargin{:});
            
            this.imageComposite_ = mlfourdfp.ImageComposite(this, 'indicesLogical', ip.Results.indicesLogical);
            this.indexOfReference = ip.Results.indexOfReference;
            this.blurArg_ = ip.Results.blurArg;  
            cd(this.sessionData.tracerLocation);           
            this.resolveTag = 'ht4r';           
        end
        
        function this = resolve(this, varargin)
            %% RESOLVE iteratively calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param destMask   "
            %  @param source     "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param t40        is the initial t4-file for the transformation:  transverse is default.
            %  @param t4         is the cumulative t4-file for the transformation.
            %  @param log        is the f.q. filename of the log file.
            %  @param useMetricGradient:  cf. ${TRANSFER}/cross_modal_intro.pdf
            %  @returns t4       is the t4-file for the transformation.
            %  @returns fqfp     is the f.q. fileprefix of the co-registered output.
            
            import mlfourdfp.*;
            ip = inputParser;
            addParameter(ip, 'dest',       {},              @iscell);
            addParameter(ip, 'source',     {},              @iscell);
            addParameter(ip, 'destMask',   'none',          @ischar);
            addParameter(ip, 'sourceMask', 'none',          @ischar);
            addParameter(ip, 'maskForImages', 'maskHetero',    @ischar);
            addParameter(ip, 'destBlur',   this.blurArg,    @isnumeric); % fwhh/mm
            addParameter(ip, 'sourceBlur', this.blurArg,    @isnumeric); % fwhh/mm
            addParameter(ip, 'indicesLogical',     this.indicesLogical,     @isnumeric); % 1 to keep; 0 to skip
            addParameter(ip, 'log',        '/dev/null',     @ischar);
            addParameter(ip, 'rnumber',    this.sessionData.rnumber, @isnumeric);
            addParameter(ip, 't40',        this.buildVisitor.transverse_t4, @(x) ischar(x) || iscell(x));
            parse(ip, varargin{:});
            ipr = ip.Results;     
            ipr = this.expandDest(ipr);
            ipr = this.expandBlurs(ipr);
            ipr.resolved = ipr.source; 
            if (~this.isfinished)
                while (ipr.rnumber <= this.NRevisions)
                    ipr = this.revise(ipr, ...
                        ipr.resolved, ...
                        this.fileprefixesRevision(ipr.dest, ipr.rnumber));
                    assert(ipr.rnumber < 10);
                end        
            end
            this.sessionData.rnumber = this.NRevisions;
            this = this.finalize(ipr);
            this.ipResults_ = ipr;
        end
        function ipr  = revise(this, ipr, source, dest)
            ipr.source = source;
            ipr.dest = dest;  
            ipr = this.copyForIterations(ipr); % copy ipr.source to ipr.dest          
            this.imageRegLog = loggerFilename(dest{1}, ...
                'func', 'HeterogeneousT4ResolveBuilder_heteroReg', ...
                'path', this.sessionData.tracerLocation);
            this.resolveLog = loggerFilename(dest{1}, ...
                'func', 'HeterogeneousT4ResolveBuilder_t4ResolveAndPaste', ...
                'path', this.sessionData.tracerLocation);
            
            ipr = this.heteroReg(ipr);
            ipr = this.t4ResolveHetero(ipr); 
            ipr = this.teardownRevision(ipr);
            ipr.rnumber = ipr.rnumber + 1;
        end       
        function ipr  = heteroReg(this, ipr)
            toRegister = ipr.dest;
            blurred = this.lazyBlurImages(ipr);
            for m = 1:length(blurred)
                for n = 1:length(blurred)
                    if (m ~= n) 
                        try                            
                            %[maskM,this] = this.maskByGauss(toRegister{m});
                            %[maskN,this] = this.maskByGauss(toRegister{n});
                            %'destMask',   maskM, ...
                            %'sourceMask', maskN, ...
                            t4 = this.buildVisitor.filenameT4(toRegister{n}, toRegister{m});
                            if (~lexist(t4))
                                this.buildVisitor.align_multiSpectral( ...
                                    'dest',       blurred{m}, ...
                                    'source',     blurred{n}, ...
                                    't4',         t4, ...
                                    't4img_4dfp', false, ...
                                    'log',        this.imageRegLog);
                            end
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            handerror(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end   
        function ipr  = heteroReg2(this, ipr)
            toRegister = ipr.dest;
            for m = 1:length(toRegister)
                for n = 1:length(toRegister)
                    if (m ~= n) 
                        try
                            [maskM,this] = this.maskByGauss(toRegister{m}, 'excludeLowerFracZ', this.lowerFracZ(m));
                            [maskN,this] = this.maskByGauss(toRegister{n}, 'excludeLowerFracZ', this.lowerFracZ(n));    
                            this.buildVisitor.align_multiSpectral( ...
                                'dest',       toRegister{m}, ...
                                'source',     toRegister{n}, ...
                                'destMask',   maskM, ...
                                'sourceMask', maskN, ...
                                'destBlur',   6.09, ...
                                'sourceBlur', 6.09, ...
                                't4',         this.buildVisitor.filenameT4(toRegister{n}, toRegister{m}), ...
                                't4img_4dfp', false, ...
                                'useMetricGradient', false, ...
                                'log',        this.imageRegLog);
                            this.buildVisitor.align_2051( ...
                                'dest',       toRegister{m}, ...
                                'source',     toRegister{n}, ...
                                'destMask',   maskM, ...
                                'sourceMask', maskN, ...
                                'destBlur',   2.03, ...
                                'sourceBlur', 2.03, ...
                                't40',        this.buildVisitor.filenameT4(toRegister{n}, toRegister{m}), ...
                                't4',         this.buildVisitor.filenameT4(toRegister{n}, toRegister{m}), ...
                                't4img_4dfp', false, ...
                                'log',        this.imageRegLog);
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            handerror(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end   
        function ipr  = heteroReg_b110(this, ipr)
            toRegister = ipr.dest;
            blurred = this.lazyBlurImages(ipr);
            for m = 1:length(blurred)
                for n = 1:length(blurred)
                    if (m ~= n) 
                        try                            
                            [maskM,this] = this.msktgen_b110_4dfp(toRegister{m});
                            [maskN,this] = this.msktgen_b110_4dfp(toRegister{n});
                            this.buildVisitor.align_multiSpectral( ...
                                'dest',       blurred{m}, ...
                                'source',     blurred{n}, ...
                                'destMask',   maskM, ...
                                'sourceMask', maskN, ...
                                't4',         this.buildVisitor.filenameT4(toRegister{n}, toRegister{m}), ...
                                't4img_4dfp', false, ...
                                'log',        this.imageRegLog);
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            handerror(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end   
        function frac = lowerFracZ(~, frame)
            if (1 == frame)
                frac = 0;
            else
                frac = 0.3;
            end
        end   
        function ipr  = t4ResolveHetero(this, ipr)
            dest = cellfun(@(x) mybasename(x), ipr.dest, 'UniformOutput', false);
            imgFns = '';
            for f = 1:length(ipr.dest)
                imgFns = [imgFns dest{f} ' ']; %#ok<AGROW>
            end
                     
            this.buildVisitor.t4_resolve(this.resolveTag, imgFns, ...
                'options', '-v -m -s', ...
                'log', this.resolveLog);
            this.t4imgDest(ipr, this.resolveTag);

            ipr.resolved = cellfun(@(x) sprintf('%s_%s', x, this.resolveTag), dest, 'UniformOutput', false); 
            movefile([this.resolveTag '.mat0'], [ipr.resolved{1} '_' datestr(now, 30) '.mat0']);
            movefile([this.resolveTag '.sub'],  [ipr.resolved{1} '_' datestr(now, 30) '.sub']);
        end
        function ipr  = teardownRevision(this, ipr)
            if (this.keepForensics); return; end
            this.teardownLogs;
            this.teardownT4s;
            
            for d = 1:length(ipr.dest)
                name = ipr.dest{d};
                name = name(1:end-1);
                for r = 1:2
                    delete(sprintf('%s%i_b*.4dfp.*', name, r));
                    delete(sprintf('%s%i_C*.4dfp.*', name, r));
                    delete(sprintf('%s%i_g*.nii.gz', name, r));
                end
            end
            delete([ipr.maskForImages '*.4dfp.*']);
            %delete([ipr.dest '.4dfp.*']);
        end
        
        %% UTILITY
        
        function ipr   = copyForIterations(this, ipr)
            if (1 == ipr.rnumber)
                for s = 1:length(ipr.source)
                    ipr.dest{s} = this.fileprefixRevision(ipr.source{s}, 1);
                    this.buildVisitor.copy_4dfp(ipr.source{s}, ipr.dest{s});
                end
                return
            end
            for d = 1:length(ipr.dest)
                this.buildVisitor.copy_4dfp( ...
                    this.fileprefixResolved( ipr.dest{d}, ipr.rnumber-1), ...
                    this.fileprefixRevision(ipr.dest{d}, ipr.rnumber));
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
        function fqfps = lazyBlurImages(this, ipr)
            %% LAZYBLURIMAGES uses specifiers in ipr; will not replace any existing indicesLogical
            %  @param ipr is a struct
            
            fqfps = cell(size(ipr.dest));
            for id = 1:length(ipr.dest)
                fqfp      = ipr.dest{id};
                fqfps{id} = this.fileprefixBlurred(fqfp, ipr.destBlur{id});
                if (~this.buildVisitor.lexist_4dfp(fqfps{b}))
                    this.buildVisitor.imgblur_4dfp(fqfp, ipr.destBlur{id});
                end
            end
        end    
        function fqfps = lazyStageImages(~, ~)
            fqfps = [];
        end    
        function dest1 = pasteImages(this, ipr, varargin)
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
    end    
    
    %% PROTECTED
    
    methods (Access = protected)
        function fqt4  = mapcycle(this, indexMin, indexMax, varargin)
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
        
        function [fp,this]   = maskByGauss(this, varargin)
            ip = inputParser;
            addRequired( ip, 'fqfp', @(x) lexist([x '.4dfp.img']));
            addParameter(ip, 'mask', ['mask_' varargin{1}], @ischar);
            addParameter(ip, 'excludeLowerFracZ', 0, @isnumeric);
            parse(ip, varargin{:});
            
            import mlfourd.* mlfourdfp.*;
            fqfp = this.buildVisitor_.gauss_4dfp(ip.Results.fqfp, 0.1);           
            ic  = ImagingContext([fqfp '.4dfp.ifh']);
            ic  = ic.numericalNiftid;
            ic  = ic.thresh(ic.dipmax/8);
            ic  = ic.binarized;
            if (ip.Results.excludeLowerFracZ > 0)
                sz = ic.size;
                interval = 1:ceil(ip.Results.excludeLowerFracZ*sz(3));
                ic.img(:,:,interval) = ...
                    zeros(sz(1), sz(2), length(interval));
            end
            fp  = ip.Results.mask;
            ic.saveas([fp '.4dfp.ifh']);
            
            this = this.pushTrash([fqfp '*.4dfp.*']);
            this = this.pushTrash([fqfp '*.nii.gz']);
        end
        function [mskt,this] = msktgen_b110_4dfp(this, fqfp)
            log = loggerFilename(mybasename(fqfp), 'func', 'HeterogenousT4ResolveBuilder_msktgen_b110_4dfp');
            atlas_ = fullfile(getenv('REFDIR'), '711-2B');
            if (~lexist(this.buildVisitor.filenameT4(fqfp, mybasename(atlas_)), 'file'))
                this.buildVisitor.lns_4dfp(atlas_);
                this.buildVisitor.align_multiSpectral( ...
                    'dest', mybasename(atlas_), ...
                    'source', fqfp, ...
                    'destBlur', this.blurArg_, ...
                    'sourceBlur', this.blurArg_, ...
                    't4img_4dfp', false, ...
                    'log', log);
            end
            this.buildVisitor.msktgen_b110_4dfp(fqfp, 'options', '-S711-2B', 'log', log);
            mskt = [fqfp '_mskt'];
            %this = this.pushTrash([mskt '*.4dfp.*']);
        end
        function [mskt,this] = msktgen2_4dfp(this, fqfp)
            log = loggerFilename(mybasename(fqfp), 'func', 'HeterogenousT4ResolveBuilder_msktgen2_4dfp');
            atlas_ = fullfile(getenv('REFDIR'), '711-2B');
            if (~lexist(this.buildVisitor.filenameT4(fqfp, mybasename(atlas_)), 'file'))
                this.buildVisitor.lns_4dfp(atlas_);
                this.buildVisitor.align_multiSpectral( ...
                    'dest', atlas_, ...
                    'source', fqfp, ...
                    'destBlur', this.blurArg_, ...
                    'sourceBlur', this.blurArg_, ...
                    't4img_4dfp', false, ...
                    'log', log);
            end
            this.buildVisitor.msktgen2_4dfp(fqfp, 'options', '-S711-2B', 'log', log);
            mskt = [fqfp '_mskt'];            
            this = this.pushTrash([mskt '*.4dfp.*']);
        end
        function this        = t4imgDest(this, ipr, tag)
            tag = mybasename(tag);
            for d = 1:length(ipr.dest)
                destFp = ipr.dest{d};
                this.buildVisitor.t4img_4dfp( ...
                    sprintf('%s_to_%s_t4', destFp, tag), ...
                    destFp, ...
                    'out', sprintf('%s_%s', destFp, tag), ...
                    'options', ['-O' ipr.dest{1}]);
            end
        end  
        function tr          = tracerACResolvedSumt(this)
            sessd = this.sessionData;
            sessd.attenuationCorrected = true;
            tr  = [this.sessionData.tracerRevision('typ', 'fqfp') 'ResolvedSumt'];
            tr0 = [           sessd.tracerRevision('typ', 'fqfp') '_on_resolved_sumt'];
            bv  = this.buildVisitor;
            if (~bv.lexist_4dfp(tr))
                assert(bv.lexist_4dfp(tr0));
                bv.lns_4dfp(tr0, tr);
            end
        end
        function tr          = tracerNACResolvedSumt(this)
            sessd = this.sessionData;
            sessd.attenuationCorrected = false;
            tr  = [this.sessionData.tracerRevision('typ', 'fqfp') 'NACResolvedSumt'];
            tr0 = [           sessd.tracerRevision('typ', 'fqfp') '_resolved_sumt'];
            bv  = this.buildVisitor;
            if (~bv.lexist_4dfp(tr))
                assert(bv.lexist_4dfp(tr0));
                bv.lns_4dfp(tr0, tr);
            end
        end
        function tr          = tracerSumtResolved(this)
            tr = this.sessionData.tracerRevision('typ', 'fqfp');
            bv = this.buildVisitor;
            if (~bv.lexist_4dfp([tr 'SumtResolved']))
                assert(bv.lexist_4dfp([tr '_sumt_on_resolved']));
                bv.lns_4dfp([tr '_sumt_on_resolved'], [tr 'SumtResolved']);
            end
            tr = [tr 'ResolvedSumt'];
        end      
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

