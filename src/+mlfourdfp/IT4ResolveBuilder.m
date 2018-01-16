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
%        finalize(this)
            %% FINALIZE
            %  @param ipr are inputParser.Results from this.resolve.
            %  @return updates this.{rnumber,product}; calls this.teardownResolve.
            %%
        logPath(this)
        onAtlasPath(this)
%        reconstituteImages(this)        
            %% RECONSTITUTEIMAGES
            %  @param ipr are inputParser.Results from this.resolve.
            %  @param optional tag is char added to filenames of reconstituted images.
            %  @return reconstituted images saved to filesystem.
            %%
        resolve(this)
            %% RESOLVE iteratively calls t4_resolve and writes a log.
            %  @param dest       is a f.q. fileprefix.
            %  @param source     "
            %  @param destMask   "
            %  @param sourceMask "
            %  @param destBlur   is the fwhm blur applied by imgblur_4dfp to dest.            
            %  @param sourceBlur is the fwhm blur applied by imgblur_4dfp to source.
            %  @param maskForImages is a f.q. fileprefix:  default maskForImages.
            %  @param indicesLogical is logical array.
            %  @param t40        is the initial t4-file for the transformation:  default this.buildVisitor.transverse_t4.
            %  @param resolveTag is char:  default this.resolveTag.
            %  @param log        is the f.q. filename of the log file:  default /dev/null.
            %  @param this.finished.isfinished == true discontinues mutations and call this.alreadyFinalized.
            %  @return intermediate filesystem objects, but N.B. actions of teardownRevision, teardownResolve, finalize.
            %%
        t4Path(this)
        teardownLogs(this)
        teardownT4s(this)
        teardownRevision(this)
        teardownResolve(this)
        
%        tf = alreadyFinalized(this)        
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

