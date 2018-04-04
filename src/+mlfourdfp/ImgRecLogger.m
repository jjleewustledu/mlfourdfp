classdef ImgRecLogger < mlpipeline.AbstractLogger
	%% IMGRECLOGGER accumulates logging strings in a CellArrayList.  It is a handle class.

	%  $Revision$
 	%  was created 20-Jan-2017 10:20:22
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64. 	
    

    properties (Constant)
        FILETYPE     = 'mlfourdfp.ImgRecLogger'
        FILETYPE_EXT = '.4dfp.img.rec'
    end
    
    properties 
        includeTimeStamp = false
    end
    
	methods 		  
 		function this = ImgRecLogger(varargin)
            this = this@mlpipeline.AbstractLogger(varargin{:});  
            this.cellArrayList_ = mlpatterns.CellArrayList; % reset  
            
            if (lexist(this.fqfilename, 'file'))
                ipr = mlfourdfp.ImgRecParser.loadx(this.fqfileprefix, this.filesuffix);
                this.cons(ipr.cellContents);
            end
        end
        
        function c = clone(this)
            %% CLONE
            %  @return c is a deep copy of a handle class
            
            c = mlfourdfp.ImgRecLogger(this);
        end
        function cons(this, preface)
            %  @param preface is char or cell array.
            %  @return internal_representation := [header preface internal_representation footer].
            
            preface = ensureCell(preface);
            newList = mlpatterns.CellArrayList; % handle
            newList.add(this.header);
            for p = 1:length(preface)
                newList.add(preface{p});
            end
            for q = 1:length(this.cellArrayList_)
                newList.add(this.cellArrayList_.get(q));
            end
            newList.add(this.footer);
            this.cellArrayList_ = clone(newList);
        end
        function consNoHeadFoot(this, preface)
            %  @param preface is char or cell array.
            %  @return internal_representation := [header preface internal_representation footer].
            
            preface = ensureCell(preface);
            newList = mlpatterns.CellArrayList; % handle
            for p = 1:length(preface)
                newList.add(preface{p});
            end
            for q = 1:length(this.cellArrayList_)
                newList.add(this.cellArrayList_.get(q));
            end
            this.cellArrayList_ = clone(newList);
        end
        function save(this)
            this = this.ensureExtension;
            mlsystem.FilesystemRegistry.cellArrayListToTextfile( ...
                this.cellArrayList_, this.fqfilename, 'w');
        end
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function fn   = defaultFqfileprefix(this)
            fn = fullfile(this.filepath, ['ImgRecLogger_' datestr(now,30)]);
        end
        function txt  = header(this)
            txt = sprintf('rec %s.4dfp.img %s %s', this.fqfileprefix, this.creationDate, this.id);
        end
        function txt  = footer(this)
            txt = sprintf('endrec %s %s', this.creationDate, this.id);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

