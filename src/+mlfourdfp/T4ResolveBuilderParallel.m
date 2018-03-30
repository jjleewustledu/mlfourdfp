classdef T4ResolveBuilderParallel < mlfourdfp.AbstractT4ResolveBuilder
	%% T4RESOLVEBUILDERPARALLEL  

	%  $Revision$
 	%  was created 10-Jan-2017 22:33:22
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties (Constant)
        CLUSTER_HOSTNAME = 'dtn01.chpc.wustl.edu'
        CLUSTER_SUBJECTS_DIR = '/scratch/jjlee/raichle/PPGdata/jjlee'
 	end

    methods (Static)
        function scp(sessFold, visit, files)
            assert(ischar(sessFold));
            if (isnumeric(visit))
                visit = sprintf('V%i', visit);
            end
            if (~iscell(files))
                files = {files};
            end
            
            import mlfourdfp.*;
            for f = 1:length(files)                
                try
                    [~,r] = mlbash(sprintf('scp -qr %s %s:%s', ...
                        fullfile(mlraichle.RaichleRegistry.instance.subjectsDir, sessFold, visit, files{f}), ...
                        T4ResolveBuilderParallel.CLUSTER_HOSTNAME, ...
                        fullfile(T4ResolveBuilderParallel.CLUSTER_SUBJECTS_DIR, sessFold, visit, '')));
                    fprintf('mlfourdfp.AbstractTracerResolveBuilder.scp:  %s\n', r);
                catch ME
                    handwarning(ME);
                end
            end
        end
        function scp_example(varargin)
            ip = inputParser;
            addOptional(ip, 'subject',  'HYGLY00', @ischar);
            addOptional(ip, 'visit',    1,         @isnumeric);
            addOptional(ip, 'patterns', {'FDG*NAC' 'mpr.4dfp.*' '*t4'}, @iscell);
            parse(ip, varargin{:});
            
            for p = 1:length(ip.Results.patterns)
                mlfourdfp.AbstractTracerResolveBuilder.scp(ip.Results.subject, ip.Results.visit, ip.Results.patterns{p});
            end
        end 
    end
    
	methods 
		  
 		function this = T4ResolveBuilderParallel(varargin)
 			%% T4RESOLVEBUILDERPARALLEL
 			%  Usage:  this = T4ResolveBuilderParallel()

 			this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});
            this.mmrBuilder_ = mlsiemens.MMRBuilder('sessionData', this.sessionData);
        end
        
        function batchClusterResolve(this)
            c = parcluster;
            c.batch(@this.batchClusterResolve__, 1, {});
        end
        function ipr   = frameRegPar(this, ipr)
            extractedFps = this.lazyStageImages(ipr);
            blurredFps   = this.lazyBlurImages(ipr);
            parfor m = this.frameRegOrdinalFrame0_:min(length(blurredFps),this.frameRegOrdinalFrameF_)
                for n = 1:length(blurredFps)
                    if (m ~= n) 
                        try
                            maskFp = this.lazyMaskForImages( ...
                                ipr.maskForImages, extractedFps{m}, extractedFps{n});
                            this.buildVisitor.align_2051( ...
                                'dest',       blurredFps{m}, ...
                                'source',     blurredFps{n}, ...
                                'destMask',   maskFp, ...
                                'sourceMask', maskFp, ...
                                't4',         this.buildVisitor.filenameT4(extractedFps{n}, extractedFps{m}), ...
                                't4img_4dfp', false, ...
                                'log',        this.imageRegLog);
                            % t4_resolve requires an idiomatic naming convention for t4 files,
                            % based on the names of frame files
                            % e. g., fdgv1r1_frame13_to_fdgv1r1_frame72_t4                            
                        catch ME
                            copyfile( ...
                                this.buildVisitor.transverse_t4, ...
                                this.buildVisitor.filenameT4(extractedFps{n}, extractedFps{m}), 'f');
                            handwarning(ME);
                        end
                    end
                end
            end            
            
            this.deleteTrash;
        end  
        function prepareClusterResolve(this)
            this.mmrBuilder_.ensureTracerLocation;
            this.pushAncillary;
            this.pushTracerNAC;
        end
        function pullTracerNAC(this, varargin)
            %% PULLTRACERNAC calls scp to pull this.CLUSTER_HOSTNAME:this.CLUSTER_SUBJECTS_DIR/<TRACER>_<VISIT>-NAC*
            %  @param visits is a cell-array defaulting to {'V1' 'V2'}
            
            ip = inputParser;
            addParameter(ip, 'visits', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            sessd = this.sessionData;
            cd(fullfile(sessd.sessionPath));
            
            for iv = 1:length(ip.Results.visits)
                cd(fullfile(sessd.sessionPath, ip.Results.visits{iv}, ...
                    sprintf('%s_%s-NAC', upper(sessd.tracer), ip.Results.visits{iv}), ''));
                listv = {'*'}; 
                for ilv = 1:length(listv)
                    try
                        mlbash(sprintf('scp -qr %s:%s .', ...
                            this.CLUSTER_HOSTNAME, ...
                            fullfile( ...
                                this.CLUSTER_SUBJECTS_DIR, sessd.sessionFolder, ip.Results.visits{iv}, ...
                                sprintf('%s_%s-NAC', upper(sessd.tracer), ip.Results.visits{iv}), listv{ilv})));                        
                    catch ME
                        handwarning(ME);
                    end
                end
            end
        end
        function pushFilesToCluster(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fqfns', @iscell);
            addRequired(ip, 'fqdn', '', @(x) lstrfind(x, 'PPGdata'));
            parse(ip, varargin{:});            
            
            fqfns = ip.Results.fqfns;
            if (~isempty(ip.Results.fqdn))
                fqdn = ip.Results.fqdn;
            else
                fqdn = fileparts(fqfns{1});
            end
            ppgDataIdx = regexp(fqdn, 'PPGdata');
            assert(~isempty(ppgDataIdx));
            fqdnCluster = sprintf('/scratch/jjlee/raichle/%s', fqdn(ppgDataIdx));
            
            mlbash(sprintf('ssh %s ''mkdir -p %s''', this.CLUSTER_HOSTNAME, fqdnCluster));
            for f = 1:length(fqfns)
                assert(lexist(fqfns{f}, 'file'));
                mlbash(sprintf('scp -qr %s %s:%s', fqfns{f}, this.CLUSTER_HOSTNAME, fqdnCluster));
            end
        end
        function pushNACUmap(this, varargin)
            ip = inputParser;
            addParameter(ip, 'visits', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            sessd = this.sessionData;
            cd(fullfile(sessd.sessionPath));            
            for iv = 1:length(ip.Results.visits)
                cd(fullfile(sessd.sessionPath, ip.Results.visits{iv}, ...
                    sprintf('%s_%s-NAC', upper(sessd.tracer), ip.Results.visits{iv}), ''));
                try
                    mlbash(sprintf('scp -qr %s %s:%s', ...
                        fullfile( ...
                            sessd.sessionPath, ip.Results.visits{iv}, ...
                            sprintf('%s_%s-NAC/%s_%s-LM-00-umap.4dfp.*', ...
                                upper(sessd.tracer), ip.Results.visits{iv}, upper(sessd.tracer), ip.Results.visits{iv})), ...
                        this.CLUSTER_HOSTNAME, ...
                        fullfile( ...
                            this.CLUSTER_SUBJECTS_DIR, sessd.sessionFolder, ip.Results.visits{iv}, '') )); 
                catch ME
                    handwarning(ME);
                end
            end
        end
        function pushTracerNAC(this, varargin)
            %% PUSHTRACERNAC calls scp to push <TRACER>_<VISIT>-NAC to this.CLUSTER_HOSTNAME:this.CLUSTER_SUBJECTS_DIR
            %  @param visits is a cell-array defaulting to {'V1' 'V2'}
            
            ip = inputParser;
            addParameter(ip, 'visits', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            sessd = this.sessionData;
            cd(fullfile(sessd.sessionPath));
            
            for iv = 1:length(ip.Results.visits)
                cd(fullfile(sessd.sessionPath, ip.Results.visits{iv}, ...
                    sprintf('%s_%s-NAC', upper(sessd.tracer), ip.Results.visits{iv}), ''));
                try
                    mlbash(sprintf('scp -qr %s %s:%s', ...
                        fullfile( ...
                            sessd.sessionPath, ip.Results.visits{iv}, ...
                            sprintf('%s_%s-NAC', upper(sessd.tracer), ip.Results.visits{iv})), ...
                        this.CLUSTER_HOSTNAME, ...
                        fullfile( ...
                            this.CLUSTER_SUBJECTS_DIR, sessd.sessionFolder, ip.Results.visits{iv}, '') )); 
                catch ME
                    handwarning(ME);
                end
            end
        end
        function pushAncillary(this, varargin)
            %% PUSHANCILLARY calls scp to push ct.4dfp.* to this.CLUSTER_HOSTNAME:this.CLUSTER_SUBJECTS_DIR
            
            ip = inputParser;
            addParameter(ip, 'visits', {'V1' 'V2'}, @iscell);
            parse(ip, varargin{:});
            
            sessd = this.sessionData;            
            cd(fullfile(sessd.sessionPath));
            lists = {'ct.4dfp.*'};
            for ils = 1:length(lists)
                try
                    mlbash(sprintf('scp -qr %s %s:%s', ...
                        lists{ils}, this.CLUSTER_HOSTNAME, fullfile(this.CLUSTER_SUBJECTS_DIR, sessd.sessionFolder, '')));
                catch ME
                    handwarning(ME);
                end
            end
            
            for iv = 1:length(ip.Results.visits)                
                cd(fullfile(sessd.sessionPath, ip.Results.visits{iv}, ''));
                listv = {'mpr.4dfp.*' '*_t4' [upper(sessd.tracer) '_*-NAC']};
                for ilv = 1:length(listv)
                    try
                        mlbash(sprintf('scp -qr %s %s:%s', ...
                            listv{ilv}, ...
                            this.CLUSTER_HOSTNAME, ...
                            fullfile(this.CLUSTER_SUBJECTS_DIR, sessd.sessionFolder, ip.Results.visits{iv}, '')));
                    catch ME
                        handwarning(ME);
                    end
                end
            end
        end   
 	end 
    
    %% PRIVATE
    
    properties (Access = private)
        mmrBuilder_
    end
    
    methods (Access = private)
        function batchClusterResolve__(this)
            this.resolveConvertedNAC;
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

