classdef FourdfpdState < mlfourdfp.ImagingState
	%% FOURDFPDSTATE has-an mlfourd.CellCompositeState 
    
	%  $Revision$
 	%  was created 25-Apr-2016 22:58:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
	properties (Dependent)
        cellComposite
        mgh
        niftic
        niftid
        numericalNiftid
 	end 

	methods %% GET
        function g = get.cellComposite(this)
            this.contexth_.changeState( ...
                mlfourd.CellCompositeState(this.concreteObj_, this.contexth_));
            g = this.contexth_.cellComposite;
        end
        function g = get.mgh(this)
            this.contexth_.changeState( ...
                mlfourd.MGHState(this.concreteObj_, this.contexth_));
            g = this.contexth_.mgh;
        end
        function g = get.niftic(this)
            this.contexth_.changeState( ...
                mlfourd.NIfTIcState(this.concreteObj_, this.contexth_));
            g = this.contexth_.niftic;
        end   
        function g = get.niftid(this)
            g = this.concreteObj_;
        end   
        function g = get.numericalNiftid(this)
            this.contexth_.changeState( ...
                mlfourd.NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            g = this.contexth_.numericalNiftid;
        end
    end

    methods
        function        addLog(this, varargin)
            this.concreteObj_.addLog(varargin{:});
        end
        function lg =   getLog(this)
            lg = this.concreteObj_.logger;
        end
        function r =    rank(this)
            r = this.concreteObj_.rank;
        end
        function this = setNoclobber(this, s)
            this.concreteObj_.noclobber = logical(s);
        end
        function tf =   sizeEq(this, varargin)
            inSize = varargin{:}.niftid.size;
            thisSize = this.concreteObj_.size;
            tf = all(thisSize(1:3) == inSize(1:3));
        end
        function tf =   sizeGt(this, varargin)
            inSize = varargin{:}.niftid.size;
            thisSize = this.concreteObj_.size;
            tf = prod(thisSize(1:3)) > prod(inSize(1:3));
        end
        function tf =   sizeLt(this, varargin)
            inSize = varargin{:}.niftid.size;
            thisSize = this.concreteObj_.size;
            tf = prod(thisSize(1:3)) < prod(inSize(1:3));
        end
        function        view(this, varargin)
            this.concreteObj_.filesuffix = '.4dfp.ifh';
            this.concreteObj_.freeview(varargin{:});
        end
        
        function this = FourdfpdState(obj, h)
            if (~isa(obj, 'mlfourd.NIfTId'))
                try
                    obj = mlfourd.NIfTId(this.dedecorateNIfTId(obj));
                catch ME
                    handexcept(ME, 'mlfourd:castingError', ...
                        'FourdfpdState.ctor does not support objects of type %s', class(obj));
                end
            end
            this.concreteObj_ = obj;
            this.contexth_ = h;
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

