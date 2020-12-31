classdef TofResolveBuilder 
	%% TOFRESOLVEBUILDER  

	%  $Revision$
 	%  was created 30-Dec-2020 00:37:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.9.0.1538559 (R2020b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        atl = '222'
 		mpr = 'T1001.4dfp.hdr'
        sessionData
        tof = 'tofp.4dfp.hdr'
    end
    
    properties (Dependent)
        mprMask
        resolveTag
        t4rb
        tofMask
        tracer
    end

	methods
        
        %% GET
        
        function g = get.mprMask(this)
            g = this.mprMask_;
        end
        function g = get.resolveTag(this)
            g = sprintf('op_%s', myfileprefix(this.tracer));
        end
        function g = get.t4rb(this)
            g = this.t4rb_;
        end
        function g = get.tofMask(this)
            g = this.tofMask_;
        end
        function g = get.tracer(this)
            g = [this.sessionData.tracerResolvedOpSubject('typ', 'fp') '_avgt.4dfp.hdr'];
        end
        
        %%
		  
 		function this = TofResolveBuilder(varargin)
 			%% TOFRESOLVEBUILDER
 			%  @param required sessionData is an mlpipelin.ISessionData.

            ip = inputParser;
            addRequired(ip, 'sessionData', @(x) isa(x, 'mlpipeline.ISessionData'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.sessionData = ipr.sessionData;
            
            this.visitor_ = mlfourdfp.FourdfpVisitor();
            %this.mprMask_ = this.buildMprMask();
            this.tofMask_ = this.buildTofMask();
            this.t4rb_ = mlfourdfp.SimpleT4ResolveBuilder( ...
                'blurArg', 4.3, ...
                'maskForImages', {'none.4dfp.hdr' this.tofMask_.filename}, ...
                'resolveTag', this.resolveTag, ...
                'theImages', {this.tracer this.tof});
        end
        function ic = buildMprMask(this)
            mprfp = myfileprefix(this.mpr);
            if ~isfile(sprintf('%s_mskt.4dfp.hdr', mprfp))
                this.visitor_.msktgenMprage(mprfp)
            end
            ic = mlfourd.ImagingContext2(sprintf('%s_mskt.4dfp.hdr', mprfp));
        end
        function ic = buildTofMask(this)
            if ~isfile([myfileprefix(this.tof) '_msk.4dfp.hdr'])                
                ic = mlfourd.ImagingContext2(this.tof);
                ic = ic.blurred(4.3);
                ic = ic.numgt(0.05);
                ic = ic / ic.dipmax;
                ic = ic * 1000;
                ic.fileprefix = [myfileprefix(this.tof) '_msk'];
                ic.save()
            else
                ic = mlfourd.ImagingContext2([myfileprefix(this.tof) '_msk.4dfp.hdr']);
            end
        end
        function t4rb = resolve(this)            
            assert(isfile(this.tracer))
            assert(isfile(this.mpr))
            assert(isfile(this.tof))
            t4rb = this.t4rb_.resolve();
        end
        function t4img_on_atl(this)
            fps = cellfun(@(x) myfileprefix(x), {this.tracer this.tof}, 'UniformOutput', false);
            fpsr = fps;
            fpsr(2:end) = cellfun(@(x) [x '_' this.resolveTag], fps(2:end), 'UniformOutput', false);
            tracerPrefix = strsplit(myfileprefix(this.tracer), '_');
            tracerPrefix = tracerPrefix{1};
            for f = 1:length(fps)
                this.visitor_.t4img_4dfp( ...
                    sprintf('%s_to_TRIO_Y_NDC_t4', tracerPrefix), ...
                    fpsr{f}, ...
                    'out', sprintf('%s_%s', fpsr{f}, this.atl), ...
                    'options', ['-O' this.atl]);
            end
        end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        mprMask_
        t4rb_   
        tofMask_
        visitor_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

