classdef CTBuilder < mlfourdfp.CompositeT4ResolveBuilder
	%% CTBUILDER  
    %  TO DO:  refactor class inheritance to composition.

	%  $Revision$
 	%  was created 14-Nov-2016 22:37:22
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	
	properties
        ct_kVp = 120
        ct_rescaleSlope = 1
        ct_rescaleIntercept = -1024
        reuseCTMasked = true
        reuseCTRescaled = true
    end
    
	methods
        function [ctOnMpr,ctToMprT4] = CT2mpr_4dfp(this, ct, varargin)
            assert(lexist(this.fourdfpImg(ct), 'file'), ...
                'mlfourdfp:RuntimeError', ...
                'CTBuilder.CT2mpr_4dfp could not find %s', this.fourdfpImg(ct)); % not necessarily in pwd            
            mpr       = this.sessionData.mpr('typ', 'fqfp');
            pth       = fileparts(mpr);
            ctToMprT4 = fullfile(pth, this.buildVisitor.filenameT4(mybasename(ct), mybasename(mpr))); 
            ctOnMpr   = fullfile(pth, [mybasename(ct) '_on_' mybasename(mpr)]);            
            if (~lexist(this.fourdfpImg(ctOnMpr)))
                ctOnMpr = this.buildVisitor.CT2mpr_4dfp(mpr, ct, ...
                    'options', ['-T' this.atlas('typ', 'fqfp')], varargin{:});
            end
            assert(lexist(ctToMprT4, 'file'));        
        end
        function ctOut = rescaleCT(this, varargin)
            ip = inputParser;
            addRequired( ip, 'ctMasked', @lexist_4dfp);
            addParameter(ip, 'ctOut', this.sessionData.ctRescaled('typ', 'fqfp'), @ischar);
            parse(ip, varargin{:});
            ctOut = ip.Results.ctOut;
            if (lexist([ctOut '.4dfp.hdr'], 'file') && this.reuseCTRescaled)
                return
            end
            
            ic = mlfourd.ImagingContext2([ip.Results.ctMasked '.4dfp.hdr']);
            ic = ic * this.ct_rescaleSlope + this.ct_rescaleIntercept;            
            ic.noclobber = false;
            ic.saveas([ctOut '.4dfp.hdr']);
        end  
        function         teardownBuildUmaps(this)
            this.teardownLogs;
            this.teardownT4s;  
            this.finished.markAsFinished( ...
                'path', this.logger.filepath, 'tag', [this.finished.tag '_' myclass(this) '_teardownBuildUmaps']);           
        end        
        function         teardownLogs(this)
            ensuredir(this.getLogPath);
            try
                movefiles('*.log', this.getLogPath); 
                movefiles('*.txt', this.getLogPath);   
                movefiles('*.lst', this.getLogPath);    
                movefiles('*.mat0', this.getLogPath);   
                movefiles('*.sub', this.getLogPath); 
            catch ME
                dispwarning(ME, 'mlfourdfp:RuntimeWarning', ...
                    'PseudoCTBuilder.teardownLogs failed to move files into %s', this.getLogPath);
            end
        end
        function         teardownT4s(this)
            if (this.keepForensics); return; end
        end
        
 		function this  = CTBuilder(varargin)
 			%% CTBuilder
 			%  Usage:  this = CTBuilder()

 			this = this@mlfourdfp.CompositeT4ResolveBuilder(varargin{:});
            this.NRevisions = 2;
            this.blurArg_ = 1.5; % per Avi, 2016oct25
        end
        
    end 
    
    %% PROTECTED
    
    methods (Access = protected)
        function tf = ctdir_exists(this)
            tf = isfolder(fullfile(this.sessionData.sessionPath, 'ct', ''));
        end
        function tf = ct4dfp_exists(this)
            tf = lexist(this.sessionData.ct('typ', '.4dfp.img'));
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

