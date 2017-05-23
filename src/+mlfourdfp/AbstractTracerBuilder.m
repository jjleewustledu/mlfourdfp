classdef (Abstract) AbstractTracerBuilder
	%% ABSTRACTTRACERBUILDER.

	%  $Revision$
 	%  was created 9-Mar-2017 15:39
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2017 John Joowon Lee. 	

    
    properties 
        sessionData
    end
    
    properties (Dependent)
        buildVisitor
        compositeResolveBuilder
        finished
        framesResolveBuilder
    end
    
    methods %% GET
        function g = get.buildVisitor(this)
            g = this.buildVisitor_;
        end
        function g = get.compositeResolveBuilder(this)
            g = this.compositeResolveBuilder_;
        end
        function g = get.finished(this)
            g = this.finished_;
        end
        function g = get.framesResolveBuilder(this)
            g = this.framesResolveBuilder_;
        end
    end
    
    methods (Static)        
        function viewStudyConverted(varargin)
            ip = inputParser;
            addParameter(ip, 'ac', false, @islogical);
            addParameter(ip, 'tracer', 'FDG', @ischar);
            parse(ip, varargin{:});
            
            fv = mlfourdfp.FourdfpVisitor;
            studyd = mlraichle.StudyData;
            cd(studyd.subjectsDir);
            subjs = mlsystem.DirTool('HYGLY*');
            for d = 1:length(subjs)
                for v = 1:2
                    try
                        sessd = mlraichle.SessionData( ...
                            'studyData', studyd, 'sessionPath', subjs.fqdns{d}, 'vnumber', v, ...
                            'tracer', ip.Results.tracer, 'ac', ip.Results.ac);
                        cd(sessd.tracerListmodeLocation);
                        if (~lexist(sessd.tracerListmodeSif('typ','fn'), 'file'))
                            fv.sif_4dfp(sessd.tracerListmodeMhdr('typ','fp'))
                        end
                    catch ME
                        handwarning(ME);
                    end
                end
            end
            for d = 1:length(subjs)
                for v = 1:2
                    try
                        sessd = mlraichle.SessionData( ...
                            'studyData', studyd, 'sessionPath', subjs.fqdns{d}, 'vnumber', v, ...
                            'tracer', ip.Results.tracer, 'ac', ip.Results.ac);
                        cd(sessd.tracerListmodeLocation);
                        ic = mlfourd.ImagingContext(sessd.tracerListmodeSif('typ','fn'));
                        ic.viewer = 'fslview';
                        ic.view;
                    catch ME
                        handwarning(ME);
                    end
                end
            end
        end
    end
    
	methods
 		function this = AbstractTracerBuilder(varargin)
 			%% AbstractTracerBuilder
 			%  Usage:  this = AbstractTracerBuilder('sessionData', theSessionData)
 			%  @param named buildVisitor
            %  @param named sessionData
 			
            %% manage parameters 
            
            import mlfourdfp.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'buildVisitor', FourdfpVisitor, @(x) isa(x, 'mlfourdfp.FourdfpVisitor'));
            addParameter(ip, 'sessionData',  [],             @(x) isa(x, 'mlpipeline.ISessionData'));
            parse(ip, varargin{:});
            
            this.buildVisitor_ = ip.Results.buildVisitor;
            this.sessionData = ip.Results.sessionData;
            %this.framesResolveBuilder_ = T4ResolveBuilder('sessionData', this.sessionData);
            %this.compositeResolveBuilder_ = CompositeT4ResolveBuilder('sessionData', this.sessionData);
        end        
 	end 
    
    %% PROTECTED
    
    methods (Access = 'protected')
        function fr = firstFortranTimeFrame(this)
            NNativeFrames = this.framesResolveBuilder.imageFrames.readLength(this.sessionData.tracerRevision('typ', 'fqfp'));
            NUmapFrames   = this.framesResolveBuilder.imageFrames.readLength(this.sessionData.tracerResolved('typ', 'fqfp'));
            fr = NNativeFrames - NUmapFrames + 1;
        end
        function fp = frameFileprefix(~, fp, fr)
            fp = sprintf('%s_frame%i', fp, fr);
        end
        function f = epochNumber(~, str)
            names = regexp(str, '\w+(-|_)(E|e)poch(?<f>\d+)', 'names');
            f = str2double(names.f);
        end
        function f = frameNumber(~, str, offset)
            names = regexp(str, '\w+(-|_)(F|f)rame(?<f>\d+)', 'names');
            f = str2double(names.f) + offset;
        end
        function pth = logPath(this)
            pth = fullfile(this.sessionData.tracerLocation, 'Log', '');
            if (~isdir(pth))
                mkdir(pth);
            end
        end
        function this = pasteFrames(this, varargin)
            
            ip = inputParser;
            addRequired(ip, 'ipr', @isstruct);
            addOptional(ip, 'tag', '', @ischar);
            parse(ip, varargin{:});
            ipr = ip.Results.ipr;
            tag = mybasename(ip.Results.tag);
            
            assert(isfield(  ipr, 'dest'));
            assert(ischar(   ipr.dest));
            assert(isfield  (ipr, 'frames'));
            assert(isnumeric(ipr.frames));
            
            pasteList = sprintf('%s_%s_paste.lst', ipr.dest, tag);
            if (lexist(pasteList)); delete(pasteList); end
            
            fid = fopen(pasteList, 'w');
            for f = 1:length(ipr.frames)
                if (ipr.frames(f))
                    fqfp = this.frameFileprefix(ipr.dest, f);
                    fprintf(fid, '%s_%s.4dfp.img\n', fqfp, tag);
                end
            end
            fclose(fid);
            this.buildVisitor.paste_4dfp(pasteList, [ipr.dest '_' tag], 'options', '-a ');
        end
    end
        
    %% PRIVATE
    
    properties (Access = private)
        buildVisitor_
        compositeResolveBuilder_
        finished_
        framesResolveBuilder_
    end    
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

