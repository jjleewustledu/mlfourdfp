classdef IT4ResolveBuilder
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
        imageRegLog
        NRevisions  
        resolveLog
        resolveTag
        theImages
    end

	methods (Abstract)
        resolve(this)
            %% RESOLVE iteratively calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param indicesLogical is logical array.
            %  @param t40        is the initial t4-file for the transformation:  default this.buildVisitor.transverse_t4.
            %  @param resolveTag is char:  default this.resolveTag.
            %  @param log        is the f.q. filename of the log file:  default /dev/null.
            %  @param this.finished.isfinished == true discontinues mutations and call this.alreadyFinalized.
            %  @return intermediate filesystem objects, but N.B. actions of teardownRevision, teardownResolve, finalize.
            
            %%
            
        finalize(this)
        teardownLogs(this)
        teardownT4s(this)
        teardownRevision(this)
        teardownResolve(this)    
        t4img_4dfp(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

