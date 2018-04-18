classdef Viewer 
	%% VIEWER  

	%  $Revision$
 	%  was created 30-Mar-2018 02:27:48 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		app = 'freeview'
    end

    methods (Static)
        function [s,r] = view(varargin)
            this = mlfourdfp.Viewer;
            targs = cell(size(varargin));
            for v = 1:length(varargin)
                targs{v} = this.interpretTarget(varargin{v});
            end
            [s,r] = mlbash(sprintf('%s %s', this.app, cell2str(targs, 'AsRow', true)));
        end
    end
    
	methods 
		  
 		function this = Viewer(varargin)
            ip = inputParser;
            addOptional(ip, 'app', 'freeview', @ischar);
            parse(ip, varargin{:});
            this.app = ip.Results.app;
 		end
    end 
    
    methods (Access = private)
        function interp = interpretTarget(this, targ)
            %  @return interp is intended for use by shell apps.  
            %  It will be char, string or cell.  Consider using cell2str(..., 'AsRow', true).
            
            fv = mlfourdfp.FourdfpVisitor;
            
            switch (class(targ))
                case 'cell'
                    interp = cellfun(@(x) this.interpretTarget(x), targ, 'UniformOutput', false);
                case {'char' 'string'}
                    
                    % 4dfp extensions
                    if (lstrfind(targ, '.4dfp.img'))
                        interp = targ;
                        return
                    end
                    if (lstrfind(targ, '.4dfp.ifh'))
                        interp = replace(targ, '.4dfp.ifh', '.4dfp.img');
                        return
                    end
                    
                    % no extension found
                    if (lexist([targ '.4dfp.img'], 'file'))
                        interp = [targ '.4dfp.img'];
                        return
                    end
                    if (lexist([targ '.nii.gz'], 'file'))
                        fv.nifti_4dfp_4([targ '.nii.gz']);
                        interp = [targ '.4dfp.img'];
                        return
                    end
                    if (lexist([targ '.nii'], 'file'))
                        fv.nifti_4dfp_4([targ '.nii']);
                        interp = [targ '.4dfp.img'];
                        return
                    end
                    
                    % NIfTI extension
                    if (lstrfind(targ, '.nii.gz') && lexist([targ '.nii.gz'], 'file'))
                        fv.nifti_4dfp_4(targ);                        
                        interp = replace(targ, '.nii.gz', '.4dfp.img');
                        return
                    end
                    if (lstrfind(targ, '.nii') && lexist([targ '.nii'], 'file'))
                        fv.nifti_4dfp_4(targ);                        
                        interp = replace(targ, '.nii', '.4dfp.img');
                        return
                    end
                    
                    % command-line options
                    interp = targ;
                otherwise
                    if (isa(targ, 'mlfourd.ImagingContext') || isa(targ, 'mlfourd.INIfTI'))
                        if (~lexist(targ.fqfilename))
                            targ.filesuffix = '.4dfp.ifh';
                            targ.save;
                        end
                        interp = this.interpretTarget(targ.fqfilename);
                        return
                    end
                    error('mlfourdfp:unsupportedSwitchcase', ...
                        'class(Viewer.interpretTarget.targ)->%s', class(targ));
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

