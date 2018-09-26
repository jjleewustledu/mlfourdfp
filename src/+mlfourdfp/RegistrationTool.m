classdef RegistrationTool < handle & mlfourd.AbstractImagingTool
	%% REGISTRATIONTOOL  

	%  $Revision$
 	%  was created 10-Aug-2018 02:39:53 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = RegistrationTool(h, varargin)
 			%% REGISTRATIONTOOL
 			%  @param .

            this = this@mlfourd.AbstractImagingTool(h, varargin{:});
            this.innerImaging_ = mlfourd.ImagingFormatContext(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

