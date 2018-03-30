classdef ImagingContext < handle
	%% IMAGINGCONTEXT provides the context for a state design pattern for imaging data.  It also 
    %  provides a facade pattern for many classes that directly represent imaging data.  It's intent  
    %  is to improve the fluent expressivity of behaviors involving imaging data.
    
	%  $Revision$
 	%  was created 25-Apr-2016 22:58:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
    
	properties (Dependent)
        filename
        filepath
        fileprefix 
        filesuffix
        fqfilename
        fqfileprefix
        fqfn
        fqfp
        noclobber
        
        fourdfpc
        fourdfpd
        
        stateTypeclass
    end
    
	methods %% GET/SET
        function f = get.filename(this)
            f = this.state_.filename;
        end
        function f = get.filepath(this)
            f = this.state_.filepath;
        end
        function f = get.fileprefix(this)
            f = this.state_.fileprefix;
        end
        function f = get.filesuffix(this)
            f = this.state_.filesuffix;
        end
        function f = get.fqfilename(this)
            f = this.state_.fqfilename;
        end
        function f = get.fqfileprefix(this)
            f = this.state_.fqfileprefix;
        end
        function f = get.fqfn(this)
            f = this.state_.fqfn;
        end
        function f = get.fqfp(this)
            f = this.state_.fqfp;
        end
        function f = get.noclobber(this)
            f = this.state_.noclobber;
        end
        
        function f = get.fourdfpc(this)
            f = this.state_.fourdfpc;
        end
        function f = get.fourdfpd(this)
            f = this.state_.fourdfpd;
        end     
        function c = get.stateTypeclass(this)
            c = class(this.state_);
        end
        
        function set.filename(this, f)
            this.state_.filename = f;
        end
        function set.filepath(this, f)
            this.state_.filepath = f;
        end        
        function set.fileprefix(this, f)
            this.state_.fileprefix = f;
        end        
        function set.filesuffix(this, f)
            this.state_.filesuffix = f;
        end        
        function set.fqfilename(this, f)
            this.state_.fqfilename = f;
        end        
        function set.fqfileprefix(this, f)
            this.state_.fqfileprefix = f;
        end        
        function set.fqfn(this, f)
            this.state_.fqfn = f;
        end        
        function set.fqfp(this, f)
            this.state_.fqfp = f;
        end
        function set.noclobber(this, f)
            this.state_.noclobber = f;
        end
    end 
    
    methods (Static)
        function this = load(obj)
            %% LOAD:  cf. ctor
            
            this = mlfourdfp.ImagingContext(obj);
        end
    end

    methods
        function      add(this, varargin)
            %% ADD
            %  @param varargin are added to a composite imaging state
            
            this.state_ = this.state_.add(varargin{:});
        end
        function      addLog(this, varargin)
            %% ADDLOG
            %  @param varargin are log entries for the imaging state
            
            this.state_.addLog(varargin{:});
        end
        function f  = char(this)
            f = char(this.state_);
        end
        function      close(this)
            this.state_.close;
        end
        function c  = createIterator(this)
            %% CREATEITERATOR
            %  @return c is an iterator for a mlpatterns.Composite instance, if any
            
            c = this.state_.createIterator;
        end
        function c  = csize(this)
            %% CSIZE
            %  @return c is the size of the imaging state when it is composite
            
            c = this.state_.csize;
        end
        function      disp(this)
            disp(this.state_);
        end
        function      ensureSaved(this)
            %% ENSURESAVED saves the imaging state as this.fqfilename on the filesystem if not already saved.
            
            if (~lexist(this.fqfilename))
                this.state_.save;
            end
        end
        function f  = find(this, varargin)
            %% FIND
            %  @param varargin are objects to find within a composite imaging state
            %  %return f is the position within the composite of the object
            
            f = this.state_.find(varargin{:});
        end
        function g  = get(this, varargin)
            %% GET
            %  @param varargin are integer locations within a composite imaging state
            %  @return g is an element of the imaging state
            
            g =  mlfourdfp.ImagingContext(this.state_.get(varargin{:}));
        end
        function g  = getLog(this)
            %% GETLOG
            %  @return g is the internal logger (handle) for the imaging state
            
            g =  this.state_.getLog;
        end
        function tf = isempty(this)
            %% ISEMPTY
            %  @return tf is boolean for state emptiness
            
            tf = this.state_.isempty;
        end
        function l  = length(this)
            %% LENGTH
            %  @return l is the length of a composite imaging state
            
            l = this.state_.length;
        end
        function      rm(this, varargin)
            %% RM
            %  @param varargin are integer locations which will be removed from the imaging state.
            
            this.state_ = this.state_.rm(varargin{:});
        end
        function      save(this)
            %% SAVE saves the imaging state as this.fqfilename on the filesystem.
            
            this.state_.save;
        end
        function      saveas(this, filename)
            %% SAVEAS saves the imaging state as this.fqfilename on the filesystem.
            %  @param filename is a string that is compatible with requirements of the filesystem;
            %  it replaces internal filename & filesystem information.

            this.state_ = this.state_.saveas(filename);
        end
        function      view(this, varargin)
            %% VIEW
            %  @param are additional filenames and other arguments to pass to the viewer.
            %  @return new window with a view of the imaging state
            %  @throws mlfourdfp:IOError
            
            this.ensureAnyFormsSaved(varargin{:});
            this.state_.view(varargin{:});
        end
        
        %% CTOR
        
        function this = ImagingContext(obj)
            %% IMAGINGCONTEXT 
            %  @param obj is imaging data:  filename, INIfTI, MGH, ImagingComponent, double, [] or 
            %  ImagingContext for copy-ctor.  
            %  @return initializes context for a state design pattern.  
            %  @throws mlfourdfp:switchCaseError, mlfourdfp:unsupportedTypeclass.
            
            import mlfourdfp.*;
            if (~exist('obj', 'var'))
                return
            end
            if (isa(obj, 'mlfourdfp.ImagingContext'))
                switch (obj.stateTypeclass)
                    case 'mlfourdfp.FourdfpcState'
                        this.state_ = FourdfpcState(obj.fourdfpc, this);
                    case 'mlfourdfp.FourdfpdState'
                        this.state_ = FourdfpdState(obj.fourdfpd, this);
                    case 'mlfourdfp.FilenameState'
                        this.state_ = FilenameState(obj.fqfilename, this);
                    otherwise
                        error('mlfourdfp:switchCaseError', ...
                              'ImagingContext.ctor.obj.stateTypeclass -> %s', obj.stateTypeclass);
                end
                return
            end
            if (isa(obj, 'mlfourdfp.Fourdfpc') || isa(obj, 'mlpatterns.CellComposite') || iscell(obj))
                this.state_ = FourdfpcState(obj, this);
                return
            end
            if (isa(obj, 'mlfourdfp.Fourdfpd'))
                this.state_ = FourdfpdState(obj, this);
                return
            end
            if (ischar(obj)) % filename need not yet exist 
                this.state_ = FilenameState(obj, this);
                return
            end
            error('mlfourdfp:unsupportedTypeclass', ...
                  'class(ctor.obj)->%s\nchar(ctor.obj)->%s', class(obj), char(obj));
        end
        function c = clone(this)
            %% CLONE simplifies calling the copy constructor.
            %  @return deep copy on new handle
            
            c = mlfourdfp.ImagingContext(this);
        end
    end
    
	methods (Hidden)
        function changeState(this, s)
            %% CHANGESTATE
            %  @param s must be an ImagingState; it replaces the current internal state.
            
            assert(isa(s, 'mlfourdfp.ImagingState'));
            this.state_ = s;
        end
    end    
    
    %% PROTECTED
    
    properties (Access = protected)
        state_
    end
    
    methods (Static, Access = protected)
        function ensureAnyFormsSaved(varargin)
            for v = 1:length(varargin)
                vobj = varargin{v};
                if (isa(vobj, 'mlfourdfp.ImagingContext'))
                    import mlfourdfp.*;
                    if (strcmp(vobj.stateTypeclass, 'mlfourdfp.FourdfpcState'))
                        ImagingContext.ensureAnyFormsSaved(vobj.fourdfpc);
                    else                        
                        ImagingContext.ensureAnyFormsSaved(vobj.fourdfpd);
                    end
                end
                if (isa(vobj, 'mlfourdfp.IFourdfpc'))
                    iter = vobj.createIterator;
                    while (iter.hasNext)
                        cached = iter.next;
                        if (~lexist(cached.fqfilename, 'file'))
                            cached.save;
                        end
                    end
                end
                if (isa(vobj, 'mlfourdfp.IFourdfpd'))
                    if (~lexist(vobj.fqfilename, 'file'))
                        vobj.save;
                    end
                end
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

