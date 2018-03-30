classdef (Abstract) IImageBuilder < mlpipeline.IImageBuilder
	%% IIMAGEBUILDER  

	%  $Revision$
 	%  was created 06-Jan-2017 18:32:42
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties (Abstract)
 	end

	methods (Abstract, Static)
        fn = fourdfpHdr
        fn = fourdfpIfh
        fn = fourdfpImg
        fn = fourdfpImgRec
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

