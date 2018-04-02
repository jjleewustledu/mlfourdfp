classdef AbstractImageComposite < mlfourdfp.IImageComposite
	%% ABSTRACTIMAGECOMPOSITE  

	%  $Revision$
 	%  was created 21-Jan-2017 17:40:16
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee.
 	

    methods (Abstract)
        len = length(this)
        ii  = nonEmptyImageIndices(this)
    end
    
    properties (Dependent)
        cImageIndices
        fortranImageIndices
        indexMin
        indexMax
        indexOfReference
        indicesLogical 
        it4ResolveBuilder
    end

    methods %% GET, SET
        function c    = get.cImageIndices(this)
            c = this.fortranImageIndices - 1;
        end
        function f    = get.fortranImageIndices(this)
            indices = 1:length(this.indicesLogical);
            f = indices(this.indicesLogical ~= 0);
        end
        function g    = get.indexMin(this)
            g = this.indexMin_;
        end
        function g    = get.indexMax(this)
            g = this.indexMax_;
        end
        function g    = get.indexOfReference(this)
            g = this.indexOfReference_;
            while (~this.indicesLogical(g))
                g = g - 1;
                if (g < 1); break; end
            end
        end
        function this = set.indexOfReference(this, s)
            assert(isnumeric(s));
            if (1 <= s && s <= length(this.indicesLogical))
                this.indexOfReference_ = s;
            end
        end
        function g    = get.indicesLogical(this)
            if (isempty(this.indicesLogical_))
                g = true;
                return
            end
            g = logical(this.indicesLogical_);
            return
            %% DEBUGGIMG
            
            neii = this.nonEmptyImageIndices;
            assert(~isempty(neii));
            g = logical(this.indicesLogical_ .* neii);
        end
        function this = set.indicesLogical(this, s)
            assert(islogical(s));
            maxIndex = this.length;
            if (isempty(s))
                this.indicesLogical_ = true(1, maxIndex);
                return
            end
            if (1 == length(s))
                if (s)
                    this.indicesLogical_ = true(1, maxIndex);
                else
                    this.indicesLogical_ = false(1, maxIndex);
                end
                return
            end
            if (length(s) > maxIndex)
                this.indicesLogical_ = s(1:maxIndex);
            end
            if (length(s) <= maxIndex)
                s_ = false(1, maxIndex);
                s_(1:length(s)) = s;
                this.indicesLogical_ = s_;
            end
        end
        function g    = get.it4ResolveBuilder(this)
            g = this.it4ResolveBuilder_;
        end
    end
    
	methods		  
 		function this = AbstractImageComposite(varargin)
 			%% ABSTRACTIMAGECOMPOSITE
 			%  Usage:  this = AbstractImageComposite()
 			
 			ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'it4rb', @(x) isa(x, 'mlfourdfp.IT4ResolveBuilder'));
            addParameter(ip, 'indexOfReference', 1, @isnumeric);
            parse(ip, varargin{:});            
            
            this.it4ResolveBuilder_ = ip.Results.it4rb;
            this.indexOfReference_ = ip.Results.indexOfReference;
 		end
 	end     
    
    %% PROTECTED
    
    properties (Access = protected)
        indexMin_
        indexMax_
        indexOfReference_
        indicesLogical_
        it4ResolveBuilder_
    end
    
    methods (Access = protected)
        function [imin,imax] = findIndexBounds(this)
            imin = 1;
            if (isempty(this.indicesLogical))
                imax = 1;
                return
            end
            if (~any(this.indicesLogical))
                imax = 1;
                return
            end
            imax = length(this.indicesLogical);
            for f = 1:length(this.indicesLogical)
                if (this.indicesLogical(f))
                    imin = f;
                    break
                end
            end
            for f = length(this.indicesLogical):-1:imin
                if (this.indicesLogical(f))
                    imax = f;
                    break
                end
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

