classdef IT4ResolveBuilder
	%% IT4RESOLVEBUILDER is an interface for t4_resolve-related builders.

	%  $Revision$
 	%  was created 11-Nov-2016 13:38:25
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.   
    
    properties (Abstract)        
        imageRegLog
        resolveLog
        resolveTag
        theImages
    end

	methods (Abstract)
        resolve(this)
            %% RESOLVE calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param indicesLogical is logical array.
            %  @param t40        is the initial t4-file for the transformation:  default this.buildVisitor.transverse_t4.
            %  @param resolveTag is char; default is this.resolveTag.
            %  @param log        is the f.q. filename of the log file; default is /dev/null.
            %  @return intermediate filesystem objects, but N.B. actions of methods teardown*(), finalize().            
            %%     
            
        finalize(this)
        teardownLogs(this)
        teardownT4s(this)
        teardownResolve(this)    
        t4img_4dfp(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

