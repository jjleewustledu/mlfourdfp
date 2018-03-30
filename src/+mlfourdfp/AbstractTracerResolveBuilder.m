classdef (Abstract) AbstractTracerResolveBuilder < mlfourdfp.T4ResolveBuilder
	%% ABSTRACTTRACERRESOLVEBUILDER.
    %  pullTracerNAC, pushTracerNAC, pushAncillary support cluster services via scp.

	%  $Revision$
 	%  was created 12-Dec-2016 19:11:06
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    
	methods
 		function this = AbstractTracerResolveBuilder(varargin)
 			%% ABSTRACTTRACERRESOLVEBUILDER
 			%  Usage:  this = AbstractTracerResolveBuilder()

 			this = this@mlfourdfp.T4ResolveBuilder(varargin{:});
            this.mmrBuilder_ = mlfourdfp.MMRResolveBuilder('sessionData', this.sessionData);
        end
        
 	end 

    %% PROTECTED
    
    properties (Access = protected)
        mmrBuilder_
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

