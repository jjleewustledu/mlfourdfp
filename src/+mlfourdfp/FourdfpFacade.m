classdef FourdfpFacade < handle
	%% FOURDFPFACADE  

	%  $Revision$
 	%  was created 10-Mar-2016 21:23:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.307022 (R2016a) Prerelease for MACI64.
 	

	properties
 		
 	end

    
    properties (Dependent)
        sessionData
        t4ResolveBuilder
    end
    
    methods %% GET
        function g = get.sessionData(this)
            g = this.sessionData_;
        end
        function g = get.t4ResolveBuilder(this)
            g = this.t4ResolveBuilder_;
        end
        function set.t4ResolveBuilder(this, s)
            assert(isa(s, 'mlfourdfp.T4ResolveBuilder'));
            this.t4ResolveBuilder_ = s;
        end
    end
    
	methods
        function this = convertMhdr(this)
            iVisits = this.sessionData_.jsreconData.createIteratorForVisits;
            while (iVisits.hasNext)
                iSeries = iVisits.createIteratorForSeries;
                while (iSeries.hasNext)
                    aSeries = iSeries.next;
                    cd(aSeries.listmodePath);
                    dt = mlsystem.DirTool('*.mhdr');
                    for d = 1:dt.length
                        this.buildVisitor_.sif_4dfp(dt.fns{d});
                    end
                end
            end
        end
        function this = convertV(this)
            iVisits = this.sessionData_.jsreconData.createIteratorForVisits;
            while (iVisits.hasNext)
                iSeries = iVisits.createIteratorForSeries;
                while (iSeries.hasNext)
                    aSeries = iSeries.next;
                    cd(aSeries.listmodePath);
                    dt = mlsystem.DirTool('*.v.hdr');
                    for d = 1:dt.length
                        this.buildVisitor_.IFhdr_to_4dfp(dt.fns{d});
                    end
                end
            end
        end
        function this = convertS(this)
            iVisits = this.sessionData_.jsreconData.createIteratorForVisits;
            while (iVisits.hasNext)
                iSeries = iVisits.createIteratorForSeries;
                while (iSeries.hasNext)
                    aSeries = iSeries.next;
                    cd(aSeries.listmodePath);
                    dt = mlsystem.DirTool('*.s.hdr');
                    for d = 1:dt.length
                        this.buildVisitor_.IFhdr_to_4dfp(dt.fns{d});
                    end
                end
            end
        end
        function this = t4ResolveSubject(this)
            if (isempty(this.t4ResolveBuilder_))
                this.t4ResolveBuilder_ = mlfourdfp.T4ResolveBuilder('sessionData', this.sessionData);
            end
            this.t4ResolveBuilder_ = this.t4ResolveBuilder_.t4ResolveSubject;
        end
 		function this = FourdfpFacade(varargin)
 			%% FOURDFPFACADE
            %  @param sessionData is an mlpipeline.SessionData specifying identifiers for the study session, including
            %  Freesurfer's recon-all results (T1.mgz is in Talairach space) and all PET data.
            %  @param t4ResolveBuilder is an mlfourdfp.T4ResolveBuilder, a builder pattern.
            %  @return this is a facade pattern for imaging alignment.

            ip = inputParser;
            addParameter(ip, 'sessionData',      [], @(x) isa(x, 'mlpipeline.SessionData'));
            addParameter(ip, 't4ResolveBuilder', [], @(x) isa(x, 'mlfourdfp.T4ResolveBuilder') || isempty(x));
            addParameter(ip, 'buildVisitor', mlfourdfp.FourdfpVisitor, ...
                                                     @(x) isa(x, 'mlfourdfp.FourdfpVisitor'));
            parse(ip, varargin{:});
            
            this.sessionData_      = ip.Results.sessionData;
            this.t4ResolveBuilder_ = ip.Results.t4ResolveBuilder;
            this.buildVisitor_     = ip.Results.buildVisitor;
        end  
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        sessionData_
        t4ResolveBuilder_
        buildVisitor_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

