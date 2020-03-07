classdef AbstractSubjectBuilder < mlfourdfp.AbstractBuilder
	%% ABSTRACTSUBJECTBUILDER  

	%  $Revision$
 	%  was created 07-May-2019 01:17:44 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = AbstractSubjectBuilder(varargin)
 			%% ABSTRACTSUBJECTBUILDER
 			%  @param census.
            %  @param subjectData is an mlpipeline.ISubjectData.

 			this = this@mlfourdfp.AbstractBuilder(varargin{:});
            ip = inputParser;
            ip.KeepUnmatched = true;            
            addParameter(ip, 'census', []);
            addParameter(ip, 'subjectData', [], @(x) isa(x, 'mlpipeline.ISubjectData'));
            parse(ip, varargin{:});            
            this.census_ = ip.Results.census;
            this.subjectData_ = ip.Results.subjectData;
 		end
 	end 
    
    %% PROTECTED
    
    properties (Access = protected)
        subjectData_
    end
    
    %% PRIVATE
    
    properties (Access = private)
        census_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

