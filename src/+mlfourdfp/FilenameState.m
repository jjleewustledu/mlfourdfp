classdef FilenameState < mlfourdfp.ImagingState
	%% FILENAMESTATE 
    %  @deprecated
    
	%  $Revision$
 	%  was created 25-Apr-2016 23:01:59
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
	properties (Dependent)
        fourdfpc
        fourdfpd
    end
    
    methods %% GET
        function f  = get.fourdfpc(this)
            this.contexth_.changeState( ...
                mlfourdfp.FourdfpcState(this.fqfilename, this.contexth_));
            f = this.contexth_.fourdfpc;
        end
        function f  = get.fourdfpd(this)
            this.contexth_.changeState( ...
                mlfourdfp.FourdfpdState(this.fqfilename, this.contexth_));
            f = this.contexth_.fourdfpd;
        end
    end
    
    methods (Static)
        function this = load(varargin)
            this = mlfourdfp.FilenameState(varargin{:});
        end
    end
    
    methods
        function        view(this)
            mlbash(sprintf('freeview %s', this.concreteObj_.fqfilename));
        end
        
        function this = FilenameState(obj, h)
            try
                obj = mlio.ConcreteIO(obj);
            catch ME
                handexcept(ME, 'mlfourdfp:castingError', ...
                    'FilenameState.load does not support objects of type %s', class(obj));
            end
            this.concreteObj_ = obj;
            this.contexth_ = h;
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

