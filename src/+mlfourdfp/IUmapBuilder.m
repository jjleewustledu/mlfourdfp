classdef (Abstract) IUmapBuilder 
	%% IUMAPBUILDER  

	%  $Revision$
 	%  was created 22-Apr-2019 22:01:39 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	properties
 		
 	end

	methods (Abstract)
        umap = buildUmap(this)
               teardownBuildUmaps(this)
        umap = umapSynth(this, varargin)
  	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

