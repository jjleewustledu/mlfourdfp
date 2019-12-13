classdef CollectionResolveBuilder < mlfourdfp.AbstractBuilder
	%% COLLECTIONRESOLVEBUILDER  

	%  $Revision$
 	%  was created 07-May-2019 14:51:48 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
    properties (Constant)
        DISABLE_FINISHFILE = true % to load t4 files into class instances
    end
    
	properties (Dependent)
 		areAligned % TODO:  -> isAligned
        compositeRB
        referenceImage
        referenceTracer   
        ReferenceTracer
        sessionData
        tracer
        t4s
        workpath
    end

    methods (Static)
        function cc = copy_4dfp_with_datetime(fn_ast)
            %% Creates copy of 4dfp, including any datetime from the original (f.-q.) filename in the copy name.
            %  @param fn_ast is filename with *s for globbing.
            %  @return ll is cell-array of sym.-links.
            
            assert(ischar(fn_ast));
            dt = mlsystem.DirTool2(fn_ast);
            cc = {};
            for f = dt.fqfns
                [~,~,x] = myfileparts(f{1});
                if strcmp(x, '.4dfp.hdr')
                    re = regexp(f{1}, '(?<TRACER>\w+)_DT(?<DT>\d+)\S+/(?<tracer>\w+)\.\S+', 'names');
                    % e.g., f{1}->$SUBJECTS_DIR/sub-S123/ses-E123/FDG_DT20190507225833.000000-Converted-AC/fdg_avgt.4dfp.img
                    if ~isempty(re)
                        [tra, lbl] = strtok(re.tracer, '_'); % e.g., tra->'fdg', lbl->'_avgt'
                        copy_prefix = sprintf('%sdt%s', tra, re.DT); % e.g., fdgdt20190507225833
                        copy = [copy_prefix lbl x];
                        cc = [cc copy]; %#ok<AGROW>
                        system(sprintf('copy_4dfp %s %s', f{1}, copy));
                    end 
                end
            end
        end
        function ll = lns_with_DateTime(fn_ast)
            %% Creates sym.-link, including any datetime from the original (f.-q.) filename in the link name.
            %  @param fn_ast is filename with *s for globbing.
            %  @return ll is cell-array of sym.-links.
            
            import mlfourdfp.CollectionResolveBuilder;
            assert(ischar(fn_ast));
            dt = mlsystem.DirTool2(fn_ast);
            ll = {};
            for f = dt.fqfns
                re = regexp(f{1}, '(?<TRACER>[A-Z]+)_DT(?<DT>\d+)\S+/(?<tracer>\w+)\.\S+', 'names');
                % e.g., f{1}->$SUBJECTS_DIR/sub-S123/ses-E123/FDG_DT20190507225833.000000-Converted-AC/fdg_avgt.4dfp.img
                if ~isempty(re)
                    [tra, lbl] = strtok(re.tracer, '_'); % e.g., tra->'fdg', lbl->'_avgt'
                    link_prefix = sprintf('%sdt%s', tra, re.DT); % e.g., fdgdt20190507225833
                    [~,~,x] = myfileparts(f{1});
                    link = [link_prefix lbl x];
                    ll = [ll link]; %#ok<AGROW>
                    if strcmp(x, '.4dfp.hdr')
                        try
                            lns_4dfp(myfileprefix(f{1}), myfileprefix(link));
                        catch ME
                            handwarning(ME);
                        end
                    end
                end                
            end
        end
        function ll = lns_with_datetime(fn_ast)
            %% Creates sym.-link, including any datetime from the original (f.-q.) filename in the link name.
            %  @param fn_ast is filename with *s for globbing.
            %  @return ll is cell-array of sym.-links.
            
            import mlfourdfp.CollectionResolveBuilder;
            assert(ischar(fn_ast));
            dt = mlsystem.DirTool2(fn_ast);
            ll = {};
            for f = dt.fqfns
                re = regexp(f{1}, '(?<tracer>[a-z]+)dt(?<dt>\d+)\S*\.4dfp\.\S+', 'names');
                % e.g., f{1}->$SUBJECTS_DIR/sub-S123/ses-E123/fdg_dt20190507225833_op_fdg_on_op_fdg_avgtr1.4dfp.img
                if ~isempty(re)
                    link_prefix = sprintf('%sdt%s', re.tracer, re.dt); % e.g., fdgdt20190507225833
                    [~,~,x] = myfileparts(f{1});
                    link = [link_prefix x];
                    ll = [ll link]; %#ok<AGROW>
                    if strcmp(x, '.4dfp.hdr')
                        try
                            lns_4dfp(myfileprefix(f{1}), myfileprefix(link));
                        catch ME
                            handwarning(ME);
                        end
                    end
                end                
            end
        end
        function prefixes = uniqueFileprefixes(imgs)
            prefixes = {};
            for im = imgs
                [~,~,x] = myfileparts(im{1});
                if strcmp(x, '.4dfp.hdr')
                    prefixes = [prefixes myfileprefix(im{1})]; %#ok<AGROW>
                end
            end
        end
    end
    
	methods 
        
        %% GET, SET
        
        function g = get.areAligned(this)
            g = this.areAligned_;
        end
        function g = get.compositeRB(this)
            g = this.compositeRB_;
        end
        function this = set.compositeRB(this, s)
            assert(isa(s, 'mlfourdfp.IT4ResolveBuilder'));
            this.compositeRB_ = s;
        end
        function g = get.referenceImage(this)
            %  @return ImagingContext.
            
            sessd = this.sessionData_;
            sessd.tracer = upper(this.referenceTracer);
            sessd.rnumber = this.rnumberOfSource_;
            g = sessd.tracerResolvedFinal('typ', 'mlfourd.ImagingContext2');
        end
        function g = get.ReferenceTracer(this)
            g = upper(this.referenceTracer_);
        end
        function this = set.ReferenceTracer(this, s)
            assert(ischar(s));
            this.referenceTracer_ = lower(s);
        end
        function g = get.referenceTracer(this)
            g = lower(this.referenceTracer_);
        end
        function this = set.referenceTracer(this, s)
            assert(ischar(s));
            this.referenceTracer_ = lower(s);
        end
        function g = get.sessionData(this)
            g = this.sessionData_;
        end
        function g = get.t4s(this)
            g = this.t4s_;
        end
        function g = get.tracer(this)
            g = this.sessionData.tracer;
        end
        function this = set.tracer(this, s)
            assert(ischar(s));
            this.sessionData_.tracer = upper(s);
        end
        function g = get.workpath(this)
            g = this.workpath_;
        end
        
        %%
        
        function this = constructFramesSubset(this, tracer, frames, varargin)
            assert(ischar(tracer));
            assert(isnumeric(frames));
            
            intermed = {};
            for d = 1:length(this.product)
                nn = this.product_{d}.nifti;
                assert(frames(end) <= size(nn, 4));
                nn.img = nn.img(:,:,:,frames);
                nn.filepath = this.workpath;
                split = strsplit(nn.fileprefix, '_');
                nn.filename = sprintf('%s_frames%ito%i.4dfp.hdr', split{1}, frames(1), frames(end));
                ic2 = mlfourd.ImagingContext2(nn);
                ic2 = ic2.timeAveraged;
                ic2.save;
                intermed{d} = ic2; %#ok<AGROW>
                % this.product{d}.fqfilename => 
                % /data/nil-bluearc/raichle/PPGdata/jjlee2/HYGLY28/Vall/fdgv1r1_op_fdgv1r1_frames1to8_sumt.4dfp.hdr

            end   
            this = this.packageProduct(intermed);
        end
        function        constructReferenceTracerToT1001T4(this)
            pwd0 = pushd(this.workpath);  
            ref  = this.referenceTracer;
            t40 = sprintf('%sr1_avgtr1_op_%sr1_avgr1_to_T1001r1_t4', ref, ref);
            if (~lexist(t40, 'file'))
                t40 = sprintf('%sr1_avgt_avgr1_to_T1001r1_t4', ref);
            end
            copyfile(t40, sprintf('%s_to_T1001_t4', ref));
            popd(pwd0);
        end
        function fp   = fileprefixStandardized(this, fp)
            toks = regexp(fp, sprintf('^(?<tracRev>\\w+)_op_\\w+_on_op_%s\\w+$', this.referenceTracer), 'names');
            if (isempty(toks))
                toks = regexp(fp, sprintf('^(?<tracRev>\\w+)_\\w*%s\\w*$', this.referenceTracer), 'names');
                assert(~isempty(toks), ...
                    'mlfourdfp:emptyRegexpTokens', 'CollectionResolveBuilder.saveStandardized');
            end
            fp = sprintf('%s_op_%s', toks.tracRev, this.referenceTracer);
        end 
        function fp   = fileprefixAvgtStandardized(this, fp)
            toks = regexp(fp, sprintf('^(?<tracRev>\\w+r\\d)_op_\\w+_on_op_%s\\w+_avgtr\\d$', this.referenceTracer), 'names');
            if (isempty(toks))
                toks = regexp(fp, sprintf('^(?<tracRev>\\w+)_\\w*%s\\w*_avgtr\\d$', this.referenceTracer), 'names');
                assert(~isempty(toks), ...
                    'mlfourdfp:emptyRegexpTokens', 'CollectionResolveBuilder.saveAvgtStandardized');
            end
            fp = sprintf('%s_avgt_op_%s', toks.tracRev, this.referenceTracer);
        end
        function this = productAverage(this, tracer)
            %  @param tracer is char.
            %  @param this.product_ is 1xN cell, N >= 1.
            %  @return this.product_ is 1x1 cell.
            
            assert(ischar(tracer));
            tracer = lower(tracer);
            avgf = this.product_{1}.fourdfp;
            if ~iscell(this.product_)
                error('mlfourdfp:TypeError', ...
                    'CollectionResolveBuilder.productAverage expected this.product_ to be a cell-array')
            end
            for p = 2:length(this.product_)
                nextf = this.product_{p}.fourdfp;
                avgf.img = avgf.img + nextf.img;
            end
            avgf.img = avgf.img / length(this.product_);
            avgf.fileprefix = sprintf('%s_avg', tracer);
            this.product_ = {mlfourd.ImagingContext2(avgf)};
            this.save;
        end 
        function this = resolve(this, fqfps_avgt, varargin)
            %  @param fqfps_avgt := cell(Nvisits, Nscans) of char fqfp.
            %  @return this.compositeRB_ := compositeT4ResolveBuilder.resolved.
            %  @return this.t4s_ := compositeT4ResolveBuilder.t4s.  See also mlfourdfp.AbstractT4ResolveBuilder.cacheT4s.
            %  @return this.product := compositeT4ResolveBuilder.product.
            %  @return this.areAligned := true.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'fqfps_avgt', @iscell);
            addParameter(ip, 'NRevisions', 1, @isnumeric);
            addParameter(ip, 'maskForImages', 'Msktgen');
            addParameter(ip, 'resolveTag', sprintf('op_%sr1', fqfps_avgt{1}), @ischar);
            addParameter(ip, 'compAlignMethod', 'align_multiSpectral', @ischar);
            addParameter(ip, 'client', 'resolve_this', @ischar);
            parse(ip, fqfps_avgt, varargin{:});
            ipr = ip.Results;
            
            fqfps_avgt = this.ensure_fqfps_avgt(fqfps_avgt);
            this.referenceTracer = fqfps_avgt{1};
            
            % trivial:  nothing to resolve
            if length(fqfps_avgt) < 2
                this.t4s_ = {fullfile(getenv('RELEASE'), 'T_t4')};
                this.product_ = mlfourd.ImagingContext2(fqfps_avgt{1});
                this.areAligned_ = true;
                return
            end
            
            % non-trivial
            pwd0 = pushd(fileparts(fqfps_avgt{1}));
            this.sessionData_.compAlignMethod = ipr.compAlignMethod;
            cRB = mlfourdfp.CompositeT4ResolveBuilder( ...
                'sessionData',   this.sessionData_, ...
                'theImages',     fqfps_avgt, ...
                'maskForImages', ipr.maskForImages, ...
                'resolveTag',    ipr.resolveTag, ...
                'NRevisions',    ipr.NRevisions, ...
                'logPath',       this.getLogPath());
            cRB.neverMarkFinished = this.DISABLE_FINISHFILE;
            cRB.ignoreFinishMark  = this.DISABLE_FINISHFILE;
            this.compositeRB_ = cRB.resolve; 
            this.t4s_ = this.compositeRB_.t4s;
            this.product_ = this.compositeRB_.product;
            this.areAligned_ = true;
            moveExisting('*.mat0', fullfile(this.workpath, 'Log', ''), 'f');
            moveExisting('*.sub',  fullfile(this.workpath, 'Log', ''), 'f');
            popd(pwd0);
        end
        function        save(this)
            for p = 1:length(this.product)
                pp = this.product{p};
                pp.filesuffix = '.4dfp.hdr';
                pp.save;
            end
        end
        function        saveAvgtStandardized(this)
            for p = 1:length(this.product)
                pp = this.product{p};
                pp = pp.timeAveraged;
                pp.fileprefix = this.fileprefixAvgtStandardized(pp.fileprefix);
                pp.filesuffix = '.4dfp.hdr';
                pp.save;                
            end          
        end
        function        saveStandardized(this)
            for p = 1:length(this.product)
                pp = this.product{p};
                pp.fileprefix = this.fileprefixStandardized(pp.fileprefix);
                pp.filesuffix = '.4dfp.hdr';
                pp.save;
            end     
        end 
        function fn   = saveThis(this, varargin) 
            ip = inputParser;
            addOptional(ip, 'client', '', @ischar);
            parse(ip, varargin{:});
            fn = sprintf('mlfourdfp_CollectionResolveBuilder_%s_this.mat', ip.Results.client);
            save(fn, 'this')
        end
        function ts   = selectT4s(this, varargin)
            %% selects t4 from this.t4s that matches the sourceTracer or destTracer.
            %  @param this.rnumber == 1 even if this.NRevisions > 1.
            %  @param named sourceTracer specifies a key for lstrfind on the source term of 
            %  this.buildVisitor_.parseFilenanmeT4(this.t4s).
            %  @param named destTracer specifies a key for the dest term, respectively.
            %  @return subset of this.t4s containing matching keys.  For no matches, return {}, never {{}}.  
            
            ip = inputParser;
            addParameter(ip, 'sourceTracer', '', @ischar);
            addParameter(ip, 'destTracer',   '', @ischar);
            parse(ip, varargin{:});
            srcTr  = lower(ip.Results.sourceTracer);
            destTr = lower(ip.Results.destTracer);
            
            % trivial case
            if (isempty(srcTr) && isempty(destTr))
                ts = {};
                return
            end   
            
            ts = cell(1,this.compositeRB.NRevisions);
            r = 1;
            ts{r} = {};
            for it = 1:length(this.t4s{r})
                [s,d] = this.buildVisitor_.parseFilenameT4(this.t4s{r}{it});
                if (lstrfind(s, srcTr) || lstrfind(d, destTr))
                    ts{r} = [ts{r} this.t4s{r}{it}];
                end
            end            
            % simplify ts = {{}} to be isomorphic to trival case
            if (1 == length(ts) && isempty(ts{1}))
                ts = {};
            end
        end
        function this = sqrt(this, varargin)
            for p = 1:length(this.product)                
                this.product_{p} = mlfourd.ImagingContext2( ...
                    this.buildVisitor_.sqrt_4dfp(this.product{p}.fqfileprefix, varargin{:}));
            end
        end
        function this = t4imgc(this, varargin)
            %% T4IMGC supports t4img_4dfp; c indicates use of cell for t4s and sources.
            %  @param required t4s     {{ char } ... }.
            %  @param required sources { ImagingContext2 ... }.
            %  @param named ref          ImagingContext2.
            %  @return this.product is { ImagingContext2 ... }.
            
            ip = inputParser;
            addRequired( ip, 't4s',                 @(x) iscell(x) && iscell(x{1}));
            addRequired( ip, 'sources',             @(x) iscell(x) && isa(x{1}, 'mlfourd.ImagingContext2'));
            addParameter(ip, 'ref', varargin{2}{1}, @(x) isa(x, 'mlfourd.ImagingContext2'));
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            if ~all(size(ipr.t4s{1}) == size(ipr.sources))
                ipr.t4s{1} = repmat(ipr.t4s{1}, size(ipr.sources));
            end
            this.product_ = cell(size(ipr.sources));
            r = 1;
            for i = 1:length(this.product_)
                try
                    fqfp = this.buildVisitor_.t4img_4dfp( ...
                        ipr.t4s{r}{i}, ipr.sources{i}.fqfileprefix, 'options', ['-O' ipr.ref.fqfileprefix]);
                    ic2 = mlfourd.ImagingContext2([fqfp '.4dfp.hdr']);
                    this.product_{i} = ic2;
                catch ME
                    dispexcept(ME)
                    % consider:   fqfp = ipr.sources{i}.fqfileprefix
                end
            end
        end    
        function this = t4imgDynamicImages(this, varargin)
            %% T4IMGDYNAMICIMAGES applies accumulated this.t4s_, typically obtained from time-averages,
            %  to dynamic sources of the time-averages specified by parameter tracer.  
            %  @param optional tracer has default := this.referenceTracer.
            %  @param {this.compositeRB_ this.t4s_} := align* method.  NRevision >= 1 is managed by this.compositeRB_. 
            %  @return this.product_ := {dynamic sources for tracer} will have a messy, but unambiguous, name.
            
            ip = inputParser;
            addOptional(ip, 'tracer', this.referenceTracer, @ischar);
            addParameter(ip, 'staging_handle', @(x) isa(x, 'function_handle')) % stageSessionScans or stageSubjectScans
            parse(ip, varargin{:});
            ipr = ip.Results;
            assert(this.areAligned);  
            assert(~isempty(this.compositeRB_)); 
            assert(~isempty(this.t4s_)); 
            
            imgs = ipr.staging_handle(ipr.tracer, '');             
            this.product_ = cell(size(imgs));
            if (length(imgs) < 2)
                this.product_{1} = mlfourd.ImagingContext2([imgs{1} '.4dfp.hdr']);
                this.product_{1}.fourdfp;
                
                toks = regexp(this.t4s_{1}{1}, '\w+_to_(?<opTag>op_[a-zA-Z0-9]+)_t4$', 'names');
                fp = this.product_{1}.fileprefix;
                if (strcmp(mybasename(this.t4s_{1}{1}), 'T_t4') || isempty(toks))
                    fp1 = sprintf('%s_op_%s', fp, this.scrubDatetime(fp));
                    this.product_{1}.fileprefix = fp1;
                    this.buildVisitor_.copyfilef_4dfp(fp, fp1);
                    return
                end
                fp1 = [fp '_' toks.opTag];
                this.product_{1}.fileprefix = fp1;
                this.buildVisitor_.copyfilef_4dfp(fp, fp1);
                return
            end
            for i = 1:length(imgs)
                try
                    this.compositeRB_ = this.compositeRB_.t4img_4dfp( ...
                        this.t4s_{1}{i}, ...
                        this.frontOfFileprefix(imgs{i}), ...
                        'ref', this.frontOfFileprefix(imgs{1})); % 'out', [this.frontOfFileprefixR1(imgs{i}) '_op_' lower(ipr.tracer)], ...
                    this.product_{i} = this.compositeRB_.product;
                catch ME
                    handwarning(ME, 'mlfourdfp:RuntimeWarning', ...
                        'CollectionResolveBuilder.t4imgDynamicImages had trouble ')
                end
            end        
        end 
        function        teardownIntermediates(~)
            %deleteExisting('*_op_*_on_op_fdg*.4dfp.*');
            deleteExisting('*_b75.4dfp*');
            deleteExisting('*_mskt.4dfp*');
        end
		  
        %% UTILITY
        
        function fqfp  = ensureLastRnumber(this, fqfp, r)
            %% ENSURELASTRNUMBER
            %  @param fqfp
            %  @param r integer
            %  @return             ${fqfp}r${r}                     if not fqfp has r[0-9]
            %  @return ${fqfp_upto_r[0-9]}r${r}${fqfp_after_r[0-9]} if fqfp has r[0-9]
            
            assert(isnumeric(r));
            
            if (iscell(fqfp))
                for f = 1:length(fqfp)
                    fqfp{f} = this.ensureLastRnumber(fqfp{f}, r);
                end
                return
            end
            
            startIdx = regexp(fqfp, 'r\d');
            if (~isempty(startIdx))
                fqfp(startIdx(end)+1) = num2str(r);
                return
            end
            fqfp = sprintf('%sr%i', fqfp, r);
        end
        function c     = extractCharFromNestedCells(this, c)
            %% EXTRACTCHARFROMNESTEDCELLS
            %  @param c is a cell or any valid argument for char().
            %  @return c is char.
            %  @throws exceptions of char.
            
            if (isempty(c))
                c = '';
                return
            end
            if (iscell(c))
                % recursion
                c = this.extractCharFromNestedCells(c{1});
                return
            end

            % basecase
            c = char(c);
        end
        function front = frontOfFileprefix(this, fps, varargin)
            %  @param fps is cell (recursive) or char (base-case).
            %  @param optional avgt is boolean.
            
            ip = inputParser;
            addOptional(ip, 'avgt', false, @islogical);
            parse(ip, varargin{:});
            
            fps = mybasename(fps);
            if (iscell(fps))
                front = cellfun(@(x) this.frontOfFileprefix(x, varargin{:}), fps, 'UniformOutput', false);
                return
            end
            
            assert(ischar(fps));
            loc = regexp(fps, '_op_\w+');
            if (isempty(loc))
                front = fps;
                return
            end
            if (~ip.Results.avgt)
                front = fps(1:loc-1);
                return
            end
            if (~strcmp(fps(loc-7:loc-2), '_avgtr'))
                front = sprintf('%s_avgt', fps(1:loc-1));
            end
        end
        function front = frontOfFileprefixR1(this, fps, varargin)
            front = this.ensureLastRnumber( ...
                this.frontOfFileprefix(fps, varargin{:}), 1);
        end
        function front = frontOfT4(~, fn)
            parts = strsplit(fn, '_');
            front = parts{1};
        end
        function s     = scrubDatetime(~, s)
            re = regexp(s, '(?<tracer>\w+)dt\d+(?<remainder>\S*)', 'names');
            assert(~isempty(re))
            s = [re.tracer re.remainder];
        end 
        
 		function this = CollectionResolveBuilder(varargin)
 			%% COLLECTIONRESOLVEBUILDER
 			%  @param sessionData.
            %  @param referenceTracer is char.
            %  @param rnumberOfSource is numeric.
 			
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.ISessionData'));
            addParameter(ip, 'referenceTracer', 'fdg', @ischar);
            addParameter(ip, 'rnumberOfSource', 2, @isnumeric);
            addParameter(ip, 'workpath', pwd, @isfolder);
            parse(ip, varargin{:});

            this.sessionData_ = ip.Results.sessionData;
            this.sessionData_.attenuationCorrected = true;
            this.rnumberOfSource_ = ip.Results.rnumberOfSource;
            this.referenceTracer_ = ip.Results.referenceTracer;
            if isempty(this.sessionData_.tracer)
                this.sessionData_.tracer = this.ReferenceTracer;
            end
            this.workpath_ = ip.Results.workpath;
 		end
 	end 
    
    %% PRIVATE
    
    properties (Access = private)
        areAligned_ = false;
        compositeRB_ 
        referenceTracer_
        rnumberOfSource_
        sessionData_
        t4s_
        workpath_
    end
    
    methods (Access = private)
        function c    = ensure_fqfps_avgt(~, c)
            %% ensures that c contains fqfileprefixes containing "_avgt" and that avgt files exist
            
            assert(~isempty(c))
            for i = 1:size(c,1)
                for j = 1:size(c,2)
                    if (isa(c{i,j}, 'mlfourd.ImagingContext2'))
                        c{i,j} = c{i,j}.fqfileprefix;
                    end
                    if ~lstrfind(c{i,j}, '_avgt')
                        if isfile([c{i,j} '_avgt.4dfp.hdr'])
                            c{i,j} = [c{i,j} '_avgt.4dfp.hdr'];
                        else
                            ic2 = mlfourd.ImagingContext2(c{i,j});
                            ic2 = ic2.timeAveraged;
                            ic2.save
                            c{i,j} = ic2.fqfileprefix;
                        end
                    end
                end
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

