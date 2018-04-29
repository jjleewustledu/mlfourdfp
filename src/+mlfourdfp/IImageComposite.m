classdef (Abstract) IImageComposite 
	%% IIMAGECOMPOSITE  

	%  $Revision$
 	%  was created 18-Jan-2017 10:39:49
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    properties (Abstract)
        cImageIndices
        fortranImageIndices
        indexMin
        indexMax
        indexOfReference
        indicesLogical
        sessionData
        theImages
    end
    
	methods (Abstract)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

