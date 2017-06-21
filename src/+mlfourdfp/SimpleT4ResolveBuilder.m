classdef SimpleT4ResolveBuilder < mlfourdfp.AbstractT4ResolveBuilder
	%% SIMPLET4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 07-Jun-2017 21:07:30 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.2.0.538062 (R2017a) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = SimpleT4ResolveBuilder(varargin)
 			%% SIMPLET4RESOLVEBUILDER
 			%  Usage:  this = SimpleT4ResolveBuilder()
 			
            this = this@mlfourdfp.AbstractT4ResolveBuilder(varargin{:}); 
            cd(this.sessionData.tracerLocation);
        end
        
        function this = resolve(this)
        end
        function this = teardownResolve(this)
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

