classdef Fourdfp < mlfourd.NIfTIdecoratorProperties
	%% FOURDFP is a NIfTIdecorator that composes an internal INIfTI object according to the decorator design pattern.
    %  It is presently a stub for future development of fourdfp tools.

	%  $Revision$
 	%  was created 24-Jan-2017 19:36:07
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	
    
    properties (Constant) 
        FILETYPE      = '4DFP'
        FILETYPE_EXT  = '.4dfp.ifh'
        FOURDFP_EXT   = '.4dfp.ifh';
        SUPPORTED_EXT = {'.4dfp.ifh' '.4dfp.hdr' '.4dfp.img' '.4dfp'};
    end
    
    methods (Static)
        function this  = load(varargin)
            %% LOAD
            %  Usage:  this = MGH.load((filename[, description]); % args passed to NIfTId
            
            import mlfourd.*;
            this = mlfourdfp.Fourdfp(NIfTId.load(varargin{:}));
        end
    end
    
	methods 
        function obj = clone(this)
            obj = mlfourdfp.Fourdfp(this.component.clone);
        end
        function       save(this)
            this.component.save;
        end
        function obj = saveas(this, fqfn)
            obj = this.clone;
            [pth,fp,x] = myfileparts(fqfn);
            if (isempty(x))
                fqfp = fullfile(pth, fp);
                obj.component_ = this.component.saveas([fqfp this.FOURDFP_EXT]);            
                obj.fvisitor_.nifti_4dfp_4(fqfp);
                obj.filesuffix = this.FOURDFP_EXT;
                deleteExisting([fqfp '.nii']);
                deleteExisting([fqfp '.nii.gz']);
                return
            end
            obj.component_ = this.component.saveas(fqfn);
        end
		  
 		function this = Fourdfp(cmp, varargin)
 			%% FOURDFP
 			%  Usage:  this = Fourdfp(NIfTIdecorator_object[, option-name, option-value, ...])
 			
            import mlfourd.*; 
            this = this@mlfourd.NIfTIdecoratorProperties(cmp, varargin{:});
            if (nargin == 1 && isa(cmp, 'mlfourdfp.Fourdfp'))
                this = this.component;
                return
            end
            this = this.append_descrip('decorated by mlfourdfp.Fourdfp');
            this.component_.filesuffix = this.FOURDFP_EXT;
            this.fvisitor_ = mlfourdfp.FourdfpVisitor;
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        fvisitor_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

