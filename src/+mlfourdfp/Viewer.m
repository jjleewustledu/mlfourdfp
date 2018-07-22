classdef Viewer < mlfourd.Viewer
	%% VIEWER prepares viewer apps with mlfourdfp.Fourdfp

	%  $Revision$
 	%  was created 30-Mar-2018 02:27:48 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.  Copyright 2018 John Joowon Lee.
    
    methods (Static)
        function [s,r] = view(varargin)
            % VIEW has this.app := 'freeview'.
            
            this = mlfourdfp.Viewer;
            [s,r] = this.aview(varargin{:});
        end
    end
    
    methods
        
 		function this = Viewer(varargin)
            this = this@mlfourd.Viewer(varargin{:});
 		end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function fqt = fqtarget(~, t)
            fqt = [myfileprefix(t) '.4dfp.img'];
        end
        function fn = constructTmp4dfp(this, fn)
            fn_ = tempFqfilename(fn);
            copyfile(fn, fn_);
            this.fv.nifti_4dfp_4(fn_); % deletes *.nii[.gz]
            fn = this.fqtarget(fn_);
        end
        function [interp,todel] = interpretTarget(this, targ)
            %  @return interp is intended for use by shell apps.  
            %  It will be char, string or cell.  Consider using cell2str(..., 'AsRow', true).
                       
            switch (class(targ))
                case 'cell'
                    [interp,todel] = cellfun(@(x) this.interpretTarget(x), targ, 'UniformOutput', false);
                case {'char' 'string'}
                    
                    % 4dfp extensions
                    if (lstrfind(targ, '.4dfp.'))
                        interp = this.fqtarget(targ);
                        todel = false;
                        return
                    end
                    
                    % no extension found
                    if (lexist([targ '.4dfp.img'], 'file'))
                        interp = this.fqtarget(targ);
                        todel = false;
                        return
                    end
                    if (lexist([targ '.nii'], 'file'))
                        interp = this.constructTmp4dfp([targ '.nii']);
                        todel = true;
                        return
                    end
                    if (lexist([targ '.nii.gz'], 'file'))
                        interp = this.constructTmp4dfp([targ '.nii.gz']);
                        todel = true;
                        return
                    end
                    
                    % NIfTI extensions .nii[.gz]
                    if (lstrfind(targ, '.nii'))
                        interp = this.constructTmp4dfp(targ);
                        todel = true;
                        return
                    end
                    
                    % command-line options
                    interp = targ;
                    todel = false;
                otherwise
                    if (isa(targ, 'mlfourd.ImagingContext'))
                        if (~lexist(targ.fqfilename))
                            targ.filesuffix = '.4dfp.ifh';
                            targ.save;
                        end
                        [interp,todel] = this.interpretTarget(targ.fqfilename);
                        return
                    end
                    error('mlfourdfp:unsupportedSwitchcase', ...
                        'class(Viewer.interpretTarget.targ)->%s', class(targ));
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

