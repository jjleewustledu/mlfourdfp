classdef DynamicUmapsBuilder < mlfourdfp.AbstractUmapResolveBuilder
	%% DynamicUmapsBuilder builds umaps for Siemens e7tools, one umap for every dynamic frame of NAC PET data. 
    %  TODO:  replace magic numbers in get.framesDynamic with invariant variables.
    %  TODO:  verify this class is obsolete and delete.
    
	%  $Revision$
 	%  was created 02-Dec-2016 13:23:54
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    properties (Dependent)
        framesDynamic
        umapSynth
    end
    
    methods %% GET
        function g = get.framesDynamic(this)
            switch (this.sessionData.tracer)
                case 'FDG'
                    g = true(1,65);
                case {'HO' 'OO'}
                    g = true(1,10);
                case 'OC'
                    g = [false(1,4) true(1,10)];
                otherwise
                    error('mlfourdfp:unsupportedSwitchCase', ...
                          'DynamicUmapsBuilder.framesDynamic.tracer->%s', this.sessionData.tracer);
            end
        end
        function g = get.umapSynth(this)
            sessdUmap = this.sessionData;
            sessdUmap.tracer = '';
            g = sessdUmap.umapSynth('typ', 'fqfp');
        end
    end

    methods (Static)
        function jobs = serialize(c)
            
            assert(isa(c, 'parallel.cluster.Generic'));
            
            hyglys = { 'HYGLY05' 'HYGLY08' 'HYGLY11' 'HYGLY24' 'HYGLY25' };
            jobs   = cell(1, length(hyglys));
            jjlee  = '/scratch/jjlee/raichle/PPGdata/jjlee';
            eSessFqdns = cellfun(@(x) fullfile(jjlee, x), hyglys, 'UniformOutput', false);
            
            for iSess = 1:length(eSessFqdns)
                try
                    jobs{iSess} = c.batch(@mlfourdfp.DynamicUmapsBuilder.serialBuildUmaps, 0, {eSessFqdns{iSess}});
                catch ME
                    handwarning(ME);
                end
            end
        end
        function parWilliam
            import mlsystem.*;
            cd(mlraichle.RaichleRegistry.instance.subjectsDir);
            dt = DirTool('*');
            dtFqdns = dt.fqdns;
            parfor idt = 1:length(dtFqdns)
                for v = 1:2
                    cd(dtFqdns{idt});
                    mlfourdfp.DynamicUmapsBuilder.serialBuildUmaps(pwd, 'iVisit', v);
                end
            end
        end
        function serialBuildUmaps(varargin)
            
            import mlfourdfp.*;
            setenv('PRINTV', '');

            ip = inputParser;
            addRequired( ip, 'sessPath', @isdir);
            addParameter(ip, 'iVisit', 1, @isnumeric);
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            cd(ip.Results.sessPath);
            DynamicUmapsBuilder.diaryv('serialBuildUmaps');
            DynamicUmapsBuilder.printv('serialBuildUmaps.ip.Results.sessPath->%s\n', ip.Results.sessPath);

            eVisit = mlsystem.DirTool(ip.Results.sessPath);
            if (mlfourdfp.T4ResolveUtilities.isVisit(eVisit.fqdns{ip.Results.iVisit}))
                eTracer = mlsystem.DirTool(eVisit.fqdns{ip.Results.iVisit});
                for iTracer = 1:length(eTracer.fqdns)
                    pth = eTracer.fqdns{iTracer};
                    DynamicUmapsBuilder.printv('serialBuildUmaps.pth:  %s\n', pth);
                    if ( mlfourdfp.T4ResolveUtilities.isTracer(pth) && ...
                         mlfourdfp.T4ResolveUtilities.isNAC(pth) && ...
                        ~mlfourdfp.T4ResolveUtilities.isEmpty(pth) && ...
                         mlfourdfp.T4ResolveUtilities.matchesTag(ip.Results.sessPath, ip.Results.tag))

                        try
                            cd(pth);
                            DynamicUmapsBuilder.printv('serialBuildUmaps:  inner try pwd->%s\n', pwd);
                            sessd = mlraichle.SessionData( ...
                                'studyData',   mlraichle.StudyData, ...
                                'sessionPath', ip.Results.sessPath, ...
                                'snumber',     mlfourdfp.T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                'tracer',      mlfourdfp.T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}));
                            this = DynamicUmapsBuilder('sessionData', sessd);
                            this = this.buildUmaps;                                     
                            save(sprintf('mlfourdfp_UmapResolveBuilder_serialBuildUmaps_this_%s.mat', mydatetimestr(now)), 'this');
                        catch ME
                            handwarning(ME);
                        end
                    end
                end
            end
        end
        function triggeringReconvertUmapsToE7Format(varargin)
            ip = inputParser;
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            
            mlfourdfp.DynamicUmapsBuilder.triggering('reconvertUmapsToE7Format', 'tag', ip.Results.tag);
        end
    end
    
	methods
        
 		function this = DynamicUmapsBuilder(varargin)
 			%% DynamicUmapsBuilder
 			%  Usage:  this = DynamicUmapsBuilder()

 			this = this@mlfourdfp.AbstractUmapResolveBuilder(varargin{:});
            assert(~isempty(this.sessionData.tracer));
            this = this.updateFinished;
        end  
        
        function [this,umaps]      = buildUmaps(this, varargin)
            [this,umaps] = this.buildDynamicUmaps(varargin{:});
        end
        function [this,umapsOpDyn] = buildDynamicUmaps(this, varargin)
            ip = inputParser;
            addParameter(ip, 'indicesLogical', this.framesDynamic, @isnumeric);
            parse(ip, varargin{:});
            this.indicesLogical = ip.Results.indicesLogical;
            pwd0 = pushd(this.sessionData.tracerNACLocation);
            
            this.keepForensics = true;
            [umapOpSumDyn,sumDyn] = this.resolveUmapOpSumDynamic(this.umapSynth);
            umapsOpDyn = this.resolveUmapsOpDynamic(umapOpSumDyn, sumDyn);
            this.reconstituteImages(struct('fqfps', umapsOpDyn, 'dest1', this.umapsOpTracer));
            this.convertUmapsToE7Format(umapsOpDyn);
            %this.keepForensics = false;
            this.teardownBuildUmaps;
            popd(pwd0);
        end
        function [umapOpSumDyn,sumDyn] = resolveUmapOpSumDynamic(this, umap)
            sumDyn = this.sumTimes(this.sessionData.tracerNACRevision('typ','fqfp'));
            this.resolveTag = ['op_' mybasename(sumDyn)];
            this = this.resolveSequence(sumDyn, this.sessionData.T1('typ','fqfp'), umap);
            umapOpSumDyn = this.product{end}.fqfp;
        end   
        function umapsOpDyn        = resolveUmapsOpDynamic(this, umapOpSumDyn, sumDyn)
            umapsOpDyn = cell(1, length(this.indicesLogical));
            dynFrames = this.lazyStageImages( ...
                struct('indicesLogical', true(1,length(this.indicesLogical)), 'dest', this.sessionData.tracerNACRevision('typ','fp')));
            [~,idxInit] = max(this.indicesLogical);
            this = this.resolveSequence(dynFrames{idxInit}, sumDyn, umapOpSumDyn);
            umapOpInit = this.product{end}.fqfp;
            framesMinusInit = this.indicesLogical;
            framesMinusInit(idxInit) = 0;
            for f = 1:length(framesMinusInit)
                if (this.indicesLogical(f))
                    this.resolveTag = ['op_' dynFrames{f}];
                    this = this.resolveSequence(dynFrames{f}, sumDyn, umapOpSumDyn);
                    umapsOpDyn{f} = this.product{end}.fqfp;
                else
                    this.resolveTag = ['op_' dynFrames{idxInit}];
                    umapsOpDyn{f} = umapOpInit;
                end
            end
        end
        
        function                     teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;
            
            %ensuredir(this.onAtlasPath);
            %movefiles(sprintf('*%s*', this.atlas('typ', 'fp')), this.onAtlasPath);
            %ensuredir(this.resolveSequenceLocation);
            %movefiles(sprintf('%s*', this.resolveSequenceTag), this.resolveSequenceLocation);            
            %delete([this.resolveSequenceTag '*_frame*4dfp*']);
            
            this.finished.markAsFinished( ...
                'path', this.logger.filepath, 'tag', [this.finished.tag '_' class(this) '_teardownBuildUmaps']); 
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
    end
    
    methods (Access = protected)
    end
    
    %% PRIVATE
    
    methods (Static, Access = private)
        function this = triggering(varargin)
            
            studyd = mlraichle.StudyData;
            
            ip = inputParser;
            addRequired( ip, 'methodName', @ischar);
            addParameter(ip, 'subjectsDir', studyd.subjectsDir, @isdir);
            addParameter(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            if (~strcmp(ip.Results.subjectsDir, studyd.subjectsDir))
                studyd.subjectsDir = ip.Results.subjectsDir;
            end
            
            import mlsystem.* mlfourdfp.* ;
            DynamicUmapsBuilder.printv('triggering.ip.Results:  %s\n', struct2str(ip.Results));            
            eSess = DirTool(ip.Results.subjectsDir);
            DynamicUmapsBuilder.printv('triggering.eSess:  %s\n', cell2str(eSess.dns));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                DynamicUmapsBuilder.printv('triggering.eVisit:  %s\n', cell2str(eVisit.dns));
                for iVisit = 1:length(eVisit.fqdns)
                    
                        eTracer = DirTool(eVisit.fqdns{iVisit});
                        DynamicUmapsBuilder.printv('triggering.eTracer:  %s\n', cell2str(eTracer.dns));
                        for iTracer = 1:length(eTracer.fqdns)

                            pth___ = eTracer.fqdns{iTracer};
                            DynamicUmapsBuilder.printv('triggering.pth:  %s\n', pth___);
                            if (T4ResolveUtilities.matchesTag(pth___, ip.Results.tag) && ...
                                T4ResolveUtilities.isVisit(pth___) && ...
                                T4ResolveUtilities.isTracer(pth___) && ...
                                T4ResolveUtilities.isNAC(pth___) && ...
                               ~T4ResolveUtilities.isEmpty(pth___))
                                try
                                    sessd = mlraichle.SessionData( ...
                                        'studyData',   studyd, ...
                                        'sessionPath', eSess.fqdns{iSess}, ...
                                        'snumber',     T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                        'tracer',      T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}));
                                    this = DynamicUmapsBuilder('sessionData', sessd);
                                    this.(ip.Results.methodName);   
                                catch ME
                                    handwarning(ME);
                                end
                            end
                        end
                end                
            end            
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

