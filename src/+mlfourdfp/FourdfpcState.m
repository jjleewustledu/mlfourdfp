classdef FourdfpcState < mlfourdfp.ImagingState
	%% FOURDFPCSTATE  
    %  @deprecated

	%  $Revision$
 	%  was created 25-Apr-2016 22:58:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
	methods 
        
        %% state changes
        
        function f = mgh(this)
            this.contexth_.changeState( ...
                mlfourd.MGHState(this.concreteObj_.get(1), this.contexth_));
            f = this.contexth_.mgh;
        end
        function g = niftid(this)
            this.contexth_.changeState( ...
                mlfourd.NIfTIdState(this.concreteObj_.get(1), this.contexth_));
            g = this.contexth_.niftid;
        end    
        function g = numericalNiftid(this)
            this.contexth_.changeState( ...
                mlfourd.NIfTIdState(this.concreteObj_.get(1), this.contexth_));
            g = this.contexth_.numericalNiftid;     
        end         
    
        %% 
        
        function this = add(this, varargin)
            this.concreteObj_ = this.concreteObj_.add(varargin{:});
        end
        function        addLog(this, varargin)
            this.concreteObj_.addLog(varargin{:});
        end
        function a    = atlas(this, varargin)
            %% ATLAS builds an atlas over the composite.
            %  @param [varargin] are any ImagingContext objects.
            %  @return a is an ImagingContext with NIfTIdState.
            
            this = this.accumulateNIfTId(varargin{:});
            
            import mlfourd.*;
            niic = this.niftic;
            a = NumericalNIfTId(niic.get(1).zeros);
            a = a.timeSummed; % reduce to 3D
            for o = 1:niic.length
                b = NumericalNIfTId(niic.get(o));
                b = b.timeSummed;
                a = a + b/dipmedian(b);
            end
            a = a.append_fileprefix('_atlas');
            a = a.append_descrip('atlas');
        end
        function c    = createIterator(this)
            c = this.concreteObj_.createIterator;
        end
        function s    = csize(this)
            s = this.concreteObj_.csize;
        end
        function f    = find(this, varargin)
            f = this.concreteObj_.find(varargin{:});
        end
        function g    = get(this, varargin)
            g = this.concreteObj_.get(varargin{:});
        end
        function tf   = isempty(this)
            tf = this.concreteObj_.isempty;
        end
        function l    = length(this)
            l = this.concreteObj_.length;
        end
        function this = rm(this, varargin)
            this.concreteObj_ = this.concreteObj_.rm(varargin{:});
        end
        function        view(this, varargin)
            niid1 = this.concreteObj_.get(1);
            fns = {};
            for f = 2:this.concreteObj_.length
                niid = this.concreteObj_.get(f);
                if (~lexist(niid.fqfilename))
                    niid.save;
                end
                fns = [fns niid.fqfileprefix '.4dfp.hdr'];
            end
            viewArgin = this.segregateForView(fns, varargin{:});
            niid1.view(viewArgin{:});
        end
        
        function this = FourdfpcState(obj, h)
            if (~isa(obj, 'mlfourd.NIfTIc'))
                try
                    obj = mlfourd.NIfTIc(obj);
                catch ME
                    handexcept(ME, 'mlfourd:castingError', ...
                        'FourdfpcState.load does not support objects of type %s', class(obj));
                end
            end
            this.concreteObj_ = obj;
            this.contexth_ = h;
        end
    end 
    
    %% PRIVATE
    
    methods (Static, Access = private)
        function viewArgs = segregateForView(fns, varargin)
            if (~isempty(varargin))
                opts = varargin;
                for o = 1:length(opts)
                    if (ischar(opts{o}))
                        if (lexist(opts{o}, 'file'))
                            fns = [fns opts{o}];
                            opts(o) = [];
                        end
                    end
                end
                viewArgs = sprintf('%s %s', ...
                    cell2str(fns, 'AsRow', true), cell2str(opts, 'AsRow', true));
                return
            end
            viewArgs = fns;
        end
    end
    
    methods (Access = private)
        function this = accumulateNIfTId(this, varargin)
            for v = 1:length(varargin)
                if (isa(varargin{v}, 'mlfourd.INIfTIc'))
                    for w = 1:length( varargin{v})
                        this.concreteObj_ = this.concreteObj_.add( ...
                            mlfourd.NIfTId(varargin{v}.get(w)));
                    end
                else
                    this.concreteObj_ = this.concreteObj_.add( ...
                        mlfourd.NIfTId(varargin{v}));
                end
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

