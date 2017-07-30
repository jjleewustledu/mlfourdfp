classdef O15CUmapResolveBuilder < mlfourdfp.AbstractUmapResolveBuilder0
	%% O15CUMAPRESOLVEBUILDER  

	%  $Revision$
 	%  was created 22-Nov-2016 16:51:36
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    methods (Static)
        function jobs = serialize(c)
            
            assert(isa(c, 'parallel.cluster.Generic'));
            
            hyglys = { 'HYGLY05' 'HYGLY08' 'HYGLY11' 'HYGLY24' 'HYGLY25' };
            jobs   = cell(1, length(hyglys));
            jjlee  = '/scratch/jjlee/raichle/PPGdata/jjlee';
            eSessFqdns = cellfun(@(x) fullfile(jjlee, x), hyglys, 'UniformOutput', false);
            
            for iSess = 1:length(eSessFqdns)
                try
                    jobs{iSess} = c.batch(@mlfourdfp.O15CUmapResolveBuilder.serialBuildUmaps2, 0, {eSessFqdns{iSess}});
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
                    mlfourdfp.O15CUmapResolveBuilder.serialBuildUmaps(pwd, 'iVisit', v);
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
            O15CUmapResolveBuilder.diaryv('serialBuildUmaps');
            O15CUmapResolveBuilder.printv('serialBuildUmaps.ip.Results.sessPath->%s\n', ip.Results.sessPath);

            eVisit = mlsystem.DirTool(ip.Results.sessPath);
            if (mlfourdfp.T4ResolveUtilities.isVisit(eVisit.fqdns{ip.Results.iVisit}))
                eTracer = mlsystem.DirTool(eVisit.fqdns{ip.Results.iVisit});
                for iTracer = 1:length(eTracer.fqdns)
                    pth = eTracer.fqdns{iTracer};
                    O15CUmapResolveBuilder.printv('serialBuildUmaps.pth:  %s\n', pth);
                    if ( mlfourdfp.T4ResolveUtilities.isTracer(pth, 'OC') && ...
                         mlfourdfp.T4ResolveUtilities.isNAC(pth) && ...
                        ~mlfourdfp.T4ResolveUtilities.isEmpty(pth) && ...
                         mlfourdfp.T4ResolveUtilities.matchesTag(ip.Results.sessPath, ip.Results.tag))

                        try
                            cd(pth);
                            O15CUmapResolveBuilder.printv('serialBuildUmaps:  inner try pwd->%s\n', pwd);
                            sessd = mlraichle.SessionData( ...
                                'studyData',   mlraichle.StudyData, ...
                                'sessionPath', ip.Results.sessPath, ...
                                'snumber',     mlfourdfp.T4ResolveUtilities.scanNumber(eTracer.dns{iTracer}), ...
                                'tracer',      mlfourdfp.T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}), ...
                                'vnumber',     mlfourdfp.T4ResolveUtilities.visitNumber(eVisit.dns{ip.Results.iVisit}));
                            this = O15CUmapResolveBuilder('sessionData', sessd);
                            this = this.buildUmaps; %#ok<NASGU>                                    
                            save(sprintf('mlfourdfp_UmapResolveBuilder_serialBuildUmaps_this_%s.mat', datestr(now, 30)), 'this');
                        catch ME
                            handwarning(ME);
                        end
                    end
                end
            end
        end
    end
    
	methods 
		  
 		function this = O15CUmapResolveBuilder(varargin)
 			%% O15CUMAPRESOLVEBUILDER
 			%  Usage:  this = O15CUmapResolveBuilder()

 			this = this@mlfourdfp.AbstractUmapResolveBuilder0(varargin{:});
            this.f18UmapResolveBuilder_ = mlfourdfp.UmapResolveBuilder(varargin{:});
            sessd = this.sessionData;
            sessd.tracer = 'FDG';
            assert(this.f18UmapResolveBuilder_.isfinished);
            cd(sessd.vLocation);
            this.finished = mlpipeline.Finished(this, 'path', this.logPath, 'tag', lower(this.sessionData.tracer));
        end
        
        function this = buildUmaps(this, varargin)
            if (~this.isfinished)
                this = this.loadSessionDataCache({'OC'});
                this = this.buildO15CUmaps(varargin{:});
            end
        end
        function this = buildO15CUmaps(this, varargin)
            
            o15Nacs = {};
            for sdc = 1:length(this.sessionDataCache)
                this.sessionData_ = this.sessionDataCache{sdc};
                o15Nacs = [o15Nacs {this.buildTracerNAC}]; %#ok<AGROW>
            end
            resolved =  this.alignO15NACs(o15Nacs);             
            [this,umaps] = this.alignUmapToNACs(o15Nacs);
            
            if (~isdir(this.sessionData.tracerNACLocation))
                mkdir(this.sessionData.tracerNACLocation);
            end
            cd(this.sessionData.tracerNACLocation);
            if (this.isfinished)
                return
            end
            this.convertUmapTo4dfp;            
            
            umap         = this.sessionData.tracerListmodeUmap('typ', 'fqfp');
            [this,umap]  = this.alignUmapToSumtResolved(umap);            
            [this,umaps] = this.alignUmapToNACFrames(umap);            
            this.reconstituteImages(struct('fqfps', umaps,   'dest1', this.umapsForNAC));
            this.reconstituteImages(struct('fqfps', o15Nacs, 'dest1', this.tracerForNAC));
            this.convertUmapsToE7Format(umaps);
            this.teardownBuildUmaps;
        end
        function tof         = tofForAlignment(this, mpr)
            [~,tof] = this.buildVisitor.align_multiSpectral( ...
                'dest', mpr, ...
                'source', this.sessionData.tof('typ', 'fqfp'));
        end 
        function [this,umap] = alignUmapToSumtResolved(this, umap)
            sumtResolved = this.sumTimes(this.sessionData.tracerNACResolved('typ', 'fqfp'));
            mpr          = this.mprOnPetSumt(this.sessionData.mpr('typ', 'fqfp'), sumtResolved);   
            umap         = this.ctOnPetSumt(umap, sumtResolved);            
            t2           = this.sessionData.t2('typ', 'fqfp');
            tof          = this.sessionData.tof('typ', 'fqfp');
            umap         = this.resolveSequence(sumtResolved, mpr, t2, tof, umap);
        end
        function               teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;
            
            %ensuredir(this.onAtlasPath);
            %movefiles(sprintf('*%s*', this.atlas('typ', 'fp')), this.onAtlasPath);
            ensuredir(this.resolveSequenceLocation);
            movefiles(sprintf('%s*', this.resolveSequenceTag), this.resolveSequenceLocation);
            
            delete([this.resolveSequenceTag '*_frame*4dfp*']);
            delete(sprintf('%s_on_*.4dfp.*', this.sessionData.ct('typ', 'fp')));  
            delete(sprintf('%s_on_*.4dfp.*', this.sessionData.ctMasked('typ', 'fp')));
            this.finished.touchFinishedMarker;
        end        
    end 
    
    %% PROTECTED
        
    methods (Access = protected)
        function dest = buildO15CFromFrames(this, sessd)
            %% BUILDCO15FROMFRAMES builds 4dfp formatted CO15 NAC images from e7 indicesLogical.
            %  See also:  mlfourdfp.FourdfpVisitor.sif_4dfp.
            
            NFrames = 10;
            
            lm      = sessd.tracerListmodeSif(         'typ', 'fqfp', 'tracer', sessd.tracer, 'snumber', sessd.snumber);
            dest    = sessd.tracerNACRevision('typ', 'fqfp', 'tracer', sessd.tracer, 'snumber', sessd.snumber);   
            destLoc = sessd.tracerNACRevision('typ', 'path', 'tracer', sessd.tracer, 'snumber', sessd.snumber);  
            
            if (this.buildVisitor.lexist_4dfp(dest))
                return
            end
            if (~this.buildVisitor.lexist_4dfp(lm))
                fprintf('mlfourdfp.O15UmapResolveBuilder.buildCO15FromFrames:  building %s\n', lm);
                pwd0 = pwd;
                cd(fileparts(lm));
                
                for fr = 0:NFrames - 1
                    this.buildVisitor.IFhdr_to_4dfp(sprintf('%s_%03i_000', lm, fr), sprintf('%s_frame%i_tosum', lm, fr+1));
                end
                ipr = struct('dest', lm, 'indicesLogical', ones(1,NFrames));
                this.pasteImageIndices(ipr, 'tosum');
                delete(sprintf('%s_frame*_tosum.4dfp.*', lm));
                sumt = this.sumTimes([lm '_tosum']);
                this.buildVisitor.move_4dfp(sumt, lm);
                delete(sprintf('%s_tosum.4dfp.*', lm))
                cd(pwd0);
            end
            if (~isdir(destLoc))
                mkdir(destLoc);
            end
            ipr = struct('dest', dest, 'source', lm, 'rnumber', 0);
            this.copySourceToDest(ipr);
        end  
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
            O15CUmapResolveBuilder.printv('triggering.ip.Results:  %s\n', struct2str(ip.Results));            
            eSess = DirTool(ip.Results.subjectsDir);
            O15CUmapResolveBuilder.printv('triggering.eSess:  %s\n', cell2str(eSess.dns));
            for iSess = 1:length(eSess.fqdns)

                eVisit = DirTool(eSess.fqdns{iSess});
                O15CUmapResolveBuilder.printv('triggering.eVisit:  %s\n', cell2str(eVisit.dns));
                for iVisit = 1:length(eVisit.fqdns)
                    
                        eTracer = DirTool(eVisit.fqdns{iVisit});
                        O15CUmapResolveBuilder.printv('triggering.eTracer:  %s\n', cell2str(eTracer.dns));
                        for iTracer = 1:length(eTracer.fqdns)

                            pth___ = eTracer.fqdns{iTracer};
                            O15CUmapResolveBuilder.printv('triggering.pth:  %s\n', pth___);
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
                                        'tracer',      T4ResolveUtilities.tracerPrefix(eTracer.dns{iTracer}), ...
                                        'vnumber',     T4ResolveUtilities.visitNumber(eVisit.dns{iVisit}));
                                    this = O15CUmapResolveBuilder('sessionData', sessd);
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

