classdef PseudoCTBuilder < mlfourdfp.CarneyUmapBuilder2
	%% PSEUDOCTBUILDER manages MR-based methods to synthesize pseudo-CT results in Hounsfield units.

	%  $Revision$
 	%  was created 23-Apr-2019 12:23:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
    properties (Dependent)
        pseudoct_nii
    end
    
	methods 
        
        %% GET
        
        function g = get.pseudoct_nii(this)
            g = this.pseudoct_nii_;
        end
        
        %%
        
        function umap = buildUmap(this, varargin)
            ip = inputParser;
            addOptional(ip, 'ct', this.pseudoct_nii_, @isfile)
            parse(ip, varargin{:})
            
            ct = myfileprefix(ip.Results.ct);
            if ~isfile([ct '.nii'])
                assert(isfile([ct '.nii.gz']), 'mlfourdfp:RuntimeError', 'PseudoCTBuilder.buildUmaps')
                gunzip(ct)
            end
            system(sprintf('%s/nifti_4dfp -4 %s.nii %s.4dfp.hdr', getenv('RELEASE'), ct, ct)) % keep *.nii
            umap = this.assembleCarneyUmap(ct);
        end
        function [ct,ctToMprT4] = CT2mpr_4dfp(this, ct, varargin)
            %% @return ct unchanged and ctToMprT4 is the indentity
                        
            mpr = this.sessionData.mpr('typ', 'fqfp');
            ctToMprT4 = fullfile(this.sessionPath, this.buildVisitor.filenameT4(mybasename(ct), mybasename(mpr))); 
            copyfile(fullfile(getenv('RELEASE'), 'T_t4'), ctToMprT4)
        end
        
 		function this = PseudoCTBuilder(varargin)
 			%% PSEUDOCTBUILDER
 			%  @param .

 			this = this@mlfourdfp.CarneyUmapBuilder2(varargin{:});
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'pseudoct_nii', 't1_pct.nii', @isfile)
            parse(ip, varargin{:})
            this.pseudoct_nii_ = ip.Results.pseudoct_nii;
            this.ct_rescaleIntercept = 0;
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        pseudoct_nii_
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

