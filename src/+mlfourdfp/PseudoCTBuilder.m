classdef PseudoCTBuilder < mlfourdfp.CTBuilder
	%% PSEUDOCTBUILDER manages MR-based methods to synthesize pseudo-CT results in Hounsfield units.

	%  $Revision$
 	%  was created 23-Apr-2019 12:23:37 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	methods 
        function [ctOnMpr,ctToMprT4] = CT2mpr_4dfp(this, ct, varargin)
            %% builds MPR and Atlas T4s de novo.
            %  @param varargin are passed to mlfourdfp.FourdfpVisitor.CT2mpr_4dfp.
            
            assert(lexist(this.fourdfpImg(ct), 'file'), ...
                'mlfourdfp:RuntimeError', ...
                'PseudoCTBuilder.CT2mpr_4dfp could not find %s', this.fourdfpImg(ct));
            mpr = this.sessionData.mpr('typ', 'fqfp');
            pth = fileparts(mpr);
            ctToMprT4 = fullfile(pth, this.buildVisitor.filenameT4(mybasename(ct), mybasename(mpr))); 
            ctOnMpr   = fullfile(pth, [mybasename(ct) '_on_' mybasename(mpr)]);   
            if (~lexist([mpr '_to_' this.atlas('typ','fp') '_t4']))
                this.buildVisitor.mpr2atl_4dfp(mpr);
            end
            if (~lexist(this.fourdfpImg(ctOnMpr)))
                ctOnMpr = this.buildVisitor.CT2mpr_4dfp(mpr, ct, ...
                    'options', ['-T' this.atlas('typ','fqfp')], varargin{:});
            end
            assert(lexist(ctToMprT4, 'file'));        
        end
        
 		function this = PseudoCTBuilder(varargin)
 			%% PSEUDOCTBUILDER
 			%  @param .

 			this = this@mlfourdfp.CTBuilder(varargin{:});
            this.ct_rescaleIntercept = -1000;
 		end
 	end 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

