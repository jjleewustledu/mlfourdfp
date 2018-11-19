classdef ImageComposite < mlfourdfp.AbstractImageComposite
	%% IMAGECOMPOSITE  

	%  $Revision$
 	%  was created 18-Jan-2017 00:20:17
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	
    properties (Dependent)
        referenceImage
        sessionData   
        sourceImage
        theImages
    end

    methods 
        
        %% GET, SET
        
        function g    = get.referenceImage(this)
            g = this.theImages_{this.indexOfReference};
        end
        function g    = get.sessionData(this)
            g = this.it4ResolveBuilder_.sessionData;
        end            
        function g    = get.sourceImage(this)
            g = this.theImages;
        end 
        function g    = get.theImages(this)
            g = this.theImages_;
        end
        
        %%
        
 		function this = ImageComposite(varargin)
 			%% IMAGECOMPOSITE
 			%  Usage:  this = ImageComposite() 			

            this = this@mlfourdfp.AbstractImageComposite(varargin{:});
 			ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'it4rb', @(x) isa(x, 'mlfourdfp.IT4ResolveBuilder'));
            addParameter(ip, 'theImages', {}, @iscell);
            addParameter(ip, 'indicesLogical', true, @islogical);
            parse(ip, varargin{:});
            
            this.theImages_ = this.it4ResolveBuilder_.ensureSafeFileprefix(ip.Results.theImages);
            this.indicesLogical = ip.Results.indicesLogical;
            [this.indexMin_,this.indexMax_] = this.findIndexBounds;
        end       
        
        function len  = length(this)
            len = length(this.theImages_);
        end
 	end 
    
    %% PROTECTED
    
    properties (Access = protected)
        theImages_
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

