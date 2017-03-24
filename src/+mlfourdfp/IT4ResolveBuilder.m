classdef IT4ResolveBuilder < mlfourdfp.IImageBuilder 
	%% IT4RESOLVEBUILDER  

	%  $Revision$
 	%  was created 11-Nov-2016 13:38:25
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

    properties (Constant)
        F_HALF_x_FWHM = 10*0.4412712 % used to differentiation args of FourdfpVisitor.{imgblur_4dfp,gauss_4dfp}
    end    
    
    properties (Abstract)
        NRevisions     
        resolveTag
        
        imageRegLog
        resolveLog
    end

	methods (Abstract)
        fileprefixBlurred(this, fqfp)
        fileprefixGaussed(this, fqfp)
        fileprefixMsk(this, fqfp)
        fileprefixMskt(this, fqfp)
        fileprefixMsktgen(this, fqfp)
        fileprefixSumt(this, fqfp)   
        logPath(this)
        onAtlasPath(this)
        resolve(this)
        t4Path(this)
        teardownLogs(this)
        teardownT4s(this)
        teardownResolve(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

