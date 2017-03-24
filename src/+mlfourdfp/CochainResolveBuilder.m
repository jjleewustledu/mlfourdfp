classdef CochainResolveBuilder < mlfourdfp.AbstractT4ResolveBuilder
	%% COCHAINRESOLVEBUILDER may form composite design patterns with T4ResolveBuilder & CompositeT4ResolveBuilder. 
    %  While dense cochains for t4_resolve are conceptually simplest, in practice, cochains with selective
    %  hierarchical structures that regularize the geometrical optimization may be more numerically 
    %  efficient/feasible.  Early implementations will use cell arrays or cell composites; the aim is to 
    %  migrate to using graphs. 
    %  See also:   GoF composite pattern.

	%  $Revision$
 	%  was created 07-Feb-2017 22:21:18
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

	properties (Dependent)
 		chain
        cochain
    end
    
    methods %% GET
        function g = get.chain(this)
            g = this.chain_;
        end
        function g = get.cochain(this)
            g = this.cochain_;
        end
    end
    
    methods (Static)
        function [lenp,len,pic] = lengthPart(f)
            pic = mlpet.PETImagingContext(f);
            len = size(pic.fourdfp, 4);
            if (len > 2^8)
                error('mlsiemens:unsupportedParamValue', 'MRBuilder.partitionDynamic.len -> %i', len);
            end
            if (len > 2^2)
                lenp = floor(sqrt(len));
            end
        end   
        function this = godo()
            import mlraichle.*;
            studyd = StudyData;
            sessd  = SessionData('studyData', studyd, 'sessionPath', ...
                     fullfile(getenv('PPG'), 'jjlee2', 'HYGLY28', ''), 'vnumber', 2);
            sessd.tracer = 'FDG';
            mmrb = mlsiemens.MMRBuilder('sessionData', sessd);
            fqfp = mmrb.sif;
            fqfp = mmrb.cropfrac(fqfp);
            cd(fileparts(fqfp));
            this = mlfourdfp.CochainResolveBuilder('sessionData', sessd, 'vendorSupport', mmrb);           
            this = this.partitionTracer;
            this = this.resolvePartitions2;
            mlfourdfp.CochainResolveBuilder.godo2;
        end
        function this = godo2()
            import mlraichle.*;
            studyd = StudyData;
            sessd  = SessionData('studyData', studyd, 'sessionPath', ...
                     fullfile(getenv('PPG'), 'jjlee2', 'HYGLY28', ''), 'vnumber', 2);
            sessd.tracer = 'FDG';
            this = mlfourdfp.CochainResolveBuilder('sessionData', sessd);
            cd(sessd.tracerLocation('typ','path'));
            this.product_ = { ...
                'fdgv1_frames1-8r2_resolved_sumt' ...
                'fdgv1_frames9-16r2_resolved_sumt' ...
                'fdgv1_frames17-24r2_resolved_sumt' ...
                'fdgv1_frames25-32r2_resolved_sumt' ...
                'fdgv1_frames33-40r2_resolved_sumt' ...
                'fdgv1_frames41-48r2_resolved_sumt' ...
                'fdgv1_frames49-56r2_resolved_sumt' ...
                'fdgv1_frames57-64r2_resolved_sumt'};
            this.product_ = cellfun(@(x) mlpet.PETImagingContext([x '.4dfp.ifh']), this.product_, 'UniformOutput', false);
            this = this.sumPartitions;   
            this = this.resolveComposite; % to umapSynth
            this = this.reorderCompositeResolved; % umap to partitions
            this = this.reorderPartitionsResolved; % partitions to frames
        end
    end

	methods
 		function this = CochainResolveBuilder(varargin)
 			%% COCHAINRESOLVEBUILDER
 			%  Usage:  this = CochainResolveBuilder()

 			this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:});            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            import mlfourdfp.* mlsiemens.*;
            addParameter(ip, 'graph', graph, @(x) isa(x, 'graph'));
            addParameter(ip, 'resolvers', {}, @iscell);
            addParameter(ip, 'vendorSupport', MMRBuilder('sessionData', this.sessionData));
            parse(ip, varargin{:});
            
            this.chain_ = graph;
            this.cochain_ = graph;
            this.resolvers_ = ResolverComposite(ip.Results.resolvers);
            this.vendorSupport_ = ip.Results.vendorSupport;
        end
         
        function this = buildNAC(this)
            this.vendorSupport_.sif;
            this.vendorSupport_.cropfrac;
            this.resolvePartitions;
            this.assemblePartitions;
        end
        function this = partitionTracer(this, varargin)
            sd = this.sessionData;            
            ip = inputParser;
            addParameter(ip, 'fqfn', sd.tracerVisit('typ', '4dfp.ifh'), @(x) lexist(x, 'file'));
            parse(ip, varargin{:}); 
            [lenp,len,pic] = this.lengthPart(ip.Results.fqfn);
            
            pic.noclobber = false;            
            pics = cell(1, lenp);
            for p = 1:lenp
                ff = pic.fourdfp;
                if (p < lenp)                    
                    ff.fileprefix = sprintf('%s_frames%i-%i', pic.fileprefix, (p-1)*lenp+1, p*lenp);
                    ff.img = ff.img(:,:,:,(p-1)*lenp+1:p*lenp);
                else
                    ff.fileprefix = sprintf('%s_frames%i-%i', pic.fileprefix, (p-1)*lenp+1, len);   
                    ff.img = ff.img(:,:,:,(p-1)*lenp+1:len);
                end
                ff.save;
                pics{p} = mlpet.PETImagingContext(ff);
            end
            this.product_ = pics;
        end
        function this = resolvePartitions2(this)
            assert(iscell(this.product_) && ...
                      isa(this.product_{1}, 'mlfourd.ImagingContext'));
            pics = this.product_;
            sd = this.sessionData;
            parfor p = 1:length(pics)
                t4rb = mlfourdfp.T4ResolveBuilder('sessionData', sd, 'theImages', pics{p}.fqfp);
                t4rb.resolve('source', pics{p}.fqfp);
                pics{p} = t4rb.product;
            end
            this.product_ = pics;
        end  
        function this = sumPartitions(this)
            assert(iscell(this.product_) && ...
                      isa(this.product_{1}, 'mlfourd.ImagingContext'));
            pics = this.product_;
            for p = 1:length(pics)
                pic = pics{p};
                pic = pic.timeSummed;
                pic.fourdfp;
                pic.save;
                pics{p} = pic;
            end
            this.product_ = pics;
        end
        function this = resolveComposite(this) 
            sd = this.sessionData;         
            assert(iscell(this.product_) && ...
                      isa(this.product_{1}, 'mlfourd.ImagingContext'));
            fv = mlfourdfp.FourdfpVisitor;
            fv.lns_4dfp(sd.T1('typ','fqfp'));
            fv.lns_4dfp(sd.umapSynth('typ','fqfp','tracer',''));
            theImages = [sd.umapSynth('typ','fp','tracer','') sd.T1('typ','fp') ...
                         cellfun(@(x) x.fileprefix, this.product_, 'UniformOutput', false)];
            ct4rb = mlfourdfp.CompositeT4ResolveBuilder('sessionData', sd, 'theImages', theImages);
            ct4rb.resolve('source', theImages);
            this.product_ = ct4rb.product;
        end
        function this = reorderCompositeResolved(this)
            prod = this.product_;            
            ct4rb = mlfourdfp.CompositeT4ResolveBuilder('sessionData', this.sessionData, 'theImages', prod);
            ipr = struct('dest', []);
            ipr.dest = [prod{3:end}];
            for p = 3:length(prod)                
                ct4rb.indexOfReference = p;
                ipr = ct4rb.resolveAndPaste(ipr);
                prod{p} = ipr.resolved;
            end
            this.product_ = prod;
        end
        function this = reorderPartitionsResolved(this)
            prod = this.product_;
            prod1 = {};
            for p = 3:length(prod)
                t4rb = mlfourdfp.T4ResolveBuilder('sessionData', this.sessionData, 'theImages', prod);
                ipr = struct('dest', []);
                ipr.dest = [prod{p}];
                lenp = this.lengthPart(ipr.dest);
                for q = 1:lenp
                    t4rb.indexOfReference = q;
                    ipr = t4rb.resolveAndPaste(ipr);
                    prod1{p,q} = ipr.resolved; %#ok<AGROW>
                end
            end            
            this.product_ = prod1;
        end
        function parts = resolvePartitions(this)
                       
            lenParts = 8;
            Nparts = floor(this.imageComposite.length/lenParts);
            parts = cell(1, Nparts);
            
            parts{1} = this.imageComposite;
            parts{1}.indicesLogical = this.indicesInterval(1, lenParts);
            parts{1}.indexOfReference = lenParts;
            this.imageComposite = parts{1};
            this = this.resolvePartition( ...
                'resolveTag', sprintf('frames%i-%i_op_frame%i', 1, lenParts, parts{1}.indexOfReference));
            
            for p = 2:Nparts
                q = (p - 1)*lenParts + 1;
                parts{p} = this.imageComposite;
                parts{p}.indicesLogical = this.indicesInterval(q, q+lenParts-1);
                parts{p}.indexOfReference = q;
                this.imageComposite = parts{p};
                if (Nparts == p)
                    this.keepForensics = false;
                end
                this = this.resolvePartition( ...
                    'resolveTag', sprintf('frames%i-%i_op_frame%i', q, q+lenParts-1, parts{p}.indexOfReference));
            end
        end
        function this  = resolvePartition(this, varargin)
            this.vendorSupport_.ensureTracerLocation;
            this.vendorSupport_.ensureTracerSymlinks;
            
            sessd  = this.sessionData;
            sessd0 = this.sessionData;
            sessd0.rnumber = sessd.rnumber - 1;
            pwd_ = pushd(sessd.tracerLocation);
            this.printv('CochainResolveBuilder.resolveRevision.pwd -> %s\n', pwd);
            this = this.resolve( ...
                'dest',      sessd.tracerRevision('typ', 'fp'), ... 
                'source',    sessd0.tracerResolved('typ', 'fp'), ...
                'indicesLogical', this.indicesLogical, ...
                varargin{:});
            popd(pwd_);
        end  
        function ii = indicesInterval(this, first, last)
            ii = false(1, this.imageComposite.length);
            for i = first:last
                ii(i) = true;
            end
        end        
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        chain_
        cochain_
        resolvers_
        vendorSupport_
    end
    
    methods (Hidden)
        function this = resolve(this)
        end
        function this = teardownResolve(this)
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

