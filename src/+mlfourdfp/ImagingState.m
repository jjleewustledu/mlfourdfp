classdef (Abstract) ImagingState < mlio.IOInterface
	%% IMAGINGSTATE is the parent class for all internal states used by ImagingContext in a state design pattern.
    %  @deprecated

	%  $Revision$
 	%  was created 25-Apr-2016 22:58:52
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	 
	properties (Abstract)
%         fourdfpc
%         fourdfpd
    end
    
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
    end
    
    methods %% GET
        function f = get.filename(this)
            f = this.concreteObj_.filename;
        end
        function f = get.filepath(this)
            f = this.concreteObj_.filepath;
        end
        function f = get.fileprefix(this)
            f = this.concreteObj_.fileprefix;
        end
        function f = get.filesuffix(this)
            f = this.concreteObj_.filesuffix;
        end
        function f = get.fqfilename(this)
            f = this.concreteObj_.fqfilename;
        end
        function f = get.fqfileprefix(this)
            f = this.concreteObj_.fqfileprefix;
        end
        function f = get.fqfn(this)
            f = this.concreteObj_.fqfn;
        end
        function f = get.fqfp(this)
            f = this.concreteObj_.fqfp;
        end
        function f = get.noclobber(this)
            f = this.concreteObj_.noclobber;
        end
        
        function this = set.filename(this, f)
            this.concreteObj_.filename = f;
        end  
        function this = set.filepath(this, f)
            this.concreteObj_.filepath = f;
        end  
        function this = set.fileprefix(this, f)
            this.concreteObj_.fileprefix = f;
        end        
        function this = set.filesuffix(this, f)
            this.concreteObj_.filesuffix = f;
        end        
        function this = set.fqfilename(this, f)
            this.concreteObj_.fqfilename = f;
        end        
        function this = set.fqfileprefix(this, f)
            this.concreteObj_.fqfileprefix = f;
        end        
        function this = set.fqfn(this, f)
            this.concreteObj_.fqfilename = f;
        end        
        function this = set.fqfp(this, f)
            this.concreteObj_.fqfileprefix = f;
        end     
        function this = set.noclobber(this, f)
            this.concreteObj_.noclobber = f;
        end
    end    
        
    methods (Static)
        function obj = dedecorateFourdfpd(obj)
            if (isa(obj, 'mlfourdfp.IFourdfpd'))
                while (isa(obj, 'mlfourdfp.IFourdfpDecorator'))
                    obj = obj.component;
                end
                return
            end
        end
    end
    
    methods
        function        add(~)
            error('mlfourd:notImplemented', 'ImagingState.add'); 
        end
        function        addLog(this, varargin)
            %% ADDLOG
            %  @param [varargin] are passed to NIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            this.contexth_.addLog(varargin{:});            
        end
        function a    = atlas(this, varargin)
            %% ATLAS
            %  @param [varargin] are passed to NIfTIcState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIcState(this.concreteObj_, this.contexth_));
            a = this.contexth_.atlas(varargin{:});
        end
        function b    = binarized(this, varargin)
            %% BINARIZED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            b = this.contexth_.binarized(varargin{:});
        end
        function b    = blurred(this, varargin)
            %% BLURRED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            b = this.contexth_.blurred(varargin{:});
        end
        function c    = char(this, varargin)
            c = this.concreteObj_.char(varargin{:});
        end
        function        close(this)
            if (~lexist(this.fqfilename, 'file'))
                this.save;
            end
            this.contexth_.changeState( ...
                mlfourd.FilenameState(this.concreteObj_, this.contexth_));
        end
        function c    = createIterator(~)
            this.contexth_.changeState( ...
                mlfourd.NIfTIcState(this.concreteObj_, this.contexth_));
            c = this.contexth_.createIterator;
        end
        function s    = csize(~)
            s = 1;
        end
        function        disp(this)
            disp(this.concreteObj_);
        end
        function d    = double(this)
            d = double(this.concreteObj_);
        end
        function f    = find(this)
            f = this.concreteObj_;
        end
        function g    = get(this)
            g = this.concreteObj_;
        end
        function g    = getLog(this, varargin)
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            g = this.contexth_.getLog(varargin{:});            
        end
        function tf   = isempty(this)
            tf = isempty(this.concreteObj_);
        end
        function l    = length(~)
            l = 1;
        end
        function r    = rank(this)
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            r = this.contexth_.rank;
        end
        function        rm(~)
            error('mlfourd:notImplemented', 'ImagingState.rm');            
        end
        function m    = masked(this, varargin)
            %% MASKED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            m = this.contexth_.masked(varargin{:});
        end
        function m    = maskedByZ(this, varargin)
            %% MASKEDBYZ
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            m = this.contexth_.maskedByZ(varargin{:});
        end
        function o    = ones(this, varargin)
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            o = this.contexth_.ones(varargin{:});
        end
        function        save(this)
            this.concreteObj_.save;
        end
        function this = saveas(this, f)
            this.concreteObj_ = this.concreteObj_.saveas(f);
        end
        function        setNoclobber(this, s)
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            this.contexth_.setNoclobber(s);
        end
        function tf   = sizeEq(this, varargin)            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            tf = this.contexth_.sizeEq(varargin{:});
        end
        function tf   = sizeGt(this, varargin)            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            tf = this.contexth_.sizeGt(varargin{:});
        end
        function tf   = sizeLt(this, varargin)            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NIfTIdState(this.concreteObj_, this.contexth_));
            tf = this.contexth_.sizeLt(varargin{:});
        end
        function s    = string(this, varargin)
            s = this.concreteObj_.string(varargin{:});
        end
        function t    = thresh(this, varargin)
            %% THRESH
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            t = this.contexth_.thresh(varargin{:});
        end
        function p    = threshp(this, varargin)
            %% THRESHP
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            p = this.contexth_.threshp(varargin{:});
        end
        function s    = timeSummed(this, varargin)
            %% TIMESUMMED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            s = this.contexth_.timeSummed(varargin{:});
        end
        function s    = timeAveraged(this, varargin)
            %% TIMEAVERAGED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            s = this.contexth_.timeAveraged(varargin{:});
        end
        function s    = timeContracted(this, varargin)
            %% TIMECONTRACTED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            s = this.contexth_.timeContracted(varargin{:});
        end
        function u    = uthresh(this, varargin)
            %% UTHRESH
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            u = this.contexth_.uthresh(varargin{:});
        end
        function p    = uthreshp(this, varargin)
            %% UTHRESHP
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            p = this.contexth_.uthreshp(varargin{:});
        end
        function s    = volumeSummed(this, varargin)
            %% VOLUMESUMMED
            %  @param [varargin] are passed to NumericalNIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            s = this.contexth_.volumeSummed(varargin{:});
        end
        function        view(this, varargin)
            %% VIEW
            %  @param [varargin] are passed to NIfTIdState after a state-change
            
            import mlfourd.*;
            this.contexth_.changeState( ...
                FilenameState(this.concreteObj_, this.contexth_));
            fqfns = cellfun(@(x) x.fqfilename, varargin);
            this.contexth_.view(fqfns{:});
        end
        function z    = zeros(this, varargin)
            import mlfourd.*;
            this.contexth_.changeState( ...
                NumericalNIfTIdState(this.concreteObj_, this.contexth_));
            z = this.contexth_.zeros(varargin{:});
        end
    end
    
    methods (Hidden)
        function this = changeState(this, s)
            this.contexth_.changeState(s);
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        contexth_
        concreteObj_
    end
    
    methods (Access = protected)
        function this = ImagingState
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

