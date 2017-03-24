classdef T4ResolveUtilities 
	%% T4RESOLVEUTILITIES  
    %  TODO:  will likely need to be moved or refactored into package mlraichle.

	%  $Revision$
 	%  was created 11-Nov-2016 13:58:00
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.1.0.441655 (R2016b) for MACI64.
 	

	properties
 		
 	end

	methods (Static)  
        function c  = cell_4dfp(fqfp)
            c = { [fqfp '.4dfp.img'] ...
                  [fqfp '.4dfp.ifh'] ...
                  [fqfp '.4dfp.hdr'] ...
                  [fqfp '.4dfp.img.rec'] };
        end
        function tf = isAC(pth)
            [~,fldr] = fileparts(pth);
            tf = (~isempty(regexp(fldr, '-AC', 'once')) || ...
                  ~isempty(regexp(fldr, '-Converted-Abs$', 'once')) || ...
                  ~isempty(regexp(fldr, '-Converted$', 'once'))) ...
                 && ...
                 ~lstrfind(lower(fldr), 'backup');
        end
        function tf = isAnyTracer(pth)
            import mlfourdfp.*;
            tf = T4ResolveUtilities.isFDG(pth) || ...
                 T4ResolveUtilities.isHO( pth) || ...
                 T4ResolveUtilities.isOC( pth) || ...
                 T4ResolveUtilities.isOO( pth);
        end
        function tf = isBackup(pth)
            [~,fldr] = fileparts(pth);
            tf = ~isempty(regexp(fldr, '-Backup$', 'once'));
        end  
        function tf = isConverted(pth)
            [~,fldr] = fileparts(pth);
            tf = ~isempty(regexp(fldr, '-Converted', 'once'));
        end  
        function tf = isConvertedAC(pth)
            [~,fldr] = fileparts(pth);
            tf = ~isempty(regexp(fldr, '-Converted-AC$', 'once'));
        end  
        function tf = isConvertedNAC(pth)
            [~,fldr] = fileparts(pth);
            tf = ~isempty(regexp(fldr, '-Converted-NAC$', 'once'));
        end    
        function tf = isEmpty(pth)
            dt = mlsystem.DirTool(pth);
            tf = isempty(dt.dns) && isempty(dt.fqfns);
        end
        function tf = isFDG(pth)
            tf = mlfourdfp.T4ResolveUtilities.isTracer(pth, 'FDG');
        end
        function tf = isHO(pth)
            tf = mlfourdfp.T4ResolveUtilities.isTracer(pth, 'HO');
        end
        function tf = isNAC(pth)
            [~,fldr] = fileparts(pth);
            tf = ~isempty(regexp(fldr, '-NAC', 'once')) && ...
                 ~lstrfind(lower(fldr), 'backup');
        end
        function tf = isNotEmpty(pth)
            tf = ~mlfourdfp.T4ResolveUtilities.isEmpty(pth);
        end
        function tf = isOC(pth)
            tf = mlfourdfp.T4ResolveUtilities.isTracer(pth, 'OC');
        end
        function tf = isOO(pth)
            tf = mlfourdfp.T4ResolveUtilities.isTracer(pth, 'OO');
        end
        function tf = isTracer(varargin)
            ip = inputParser;
            addRequired(ip, 'path', @isdir);
            addOptional(ip, 'tracers', 'FDG', @(x) ischar(x) || iscell(x));
            parse(ip, varargin{:});
            
            [~,fldr] = fileparts(ip.Results.path);
            tf = lstrfind(fldr, ip.Results.tracers);
        end
        function tf = isVisit(fldr)
            tf = ~isempty(regexp(fldr, 'V[0-9]', 'once'));
        end
        function tf = matchesAC(ac, pth)
            assert(islogical(ac));
            import mlfourdfp.*;
            if (ac)
                tf =  T4ResolveUtilities.isAC(pth) && ~T4ResolveUtilities.isNAC(pth);
            else
                tf = ~T4ResolveUtilities.isAC(pth) &&  T4ResolveUtilities.isNAC(pth);
            end
        end
        function tf = matchesTag(sess, tag)
            assert(ischar(sess));
            assert(ischar(tag));
            
            if (isempty(tag))
                tf = true;
                return
            end
            tf = lstrfind(sess, tag);
        end
        function tf = pathConditions(varargin)
            ip = inputParser;
            addRequired(ip, 'pth', @isdir);
            addOptional(ip, 'conditions', @iscell);
            parse(ip, varargin{:});
            
            tf = true;
            for c = 1:length(ip.Results.conditions)
                tf = tf && mlfourdfp.T4ResolveUtilities.(ip.Results.conditions{c})(ip.Results.pth);
            end
        end
        function [s,m,mid] = safeMkdir(d)
            if (isdir(d))
                movefile(d, [d '-Backup']); 
            end
            [s,m,mid] = mkdir(d);
        end
        function [s,m,mid] = safeMovefile(f, g)
            if (isdir(g))
                movefile(g, [g '-Backup']); 
            end
            if (2 == exist(g, 'file'))
                [h1,h2,h3] = myfileparts(g);
                movefile(g, fullfile(h1, [h2 '-backup' h3]));
            end
            [s,m,mid] = movefile(f, g);
        end
        function s  = scanNumber(fldr)
            idx = regexp(upper(fldr), 'FDG|HO|OO|OC', 'end');
            s = str2double(fldr(idx+1));
        end
        function t  = tracerPrefix(varargin)
            ip = inputParser;
            addRequired( ip, 'folder', @ischar);
            addOptional( ip, 'tracers', {'FDG' 'HO' 'OO' 'OC'}, @(x) ischar(x) || iscell(x));
            addParameter(ip, 'prefixRegexp', 'FDG|HO|OO|OC', @ischar);
            parse(ip, varargin{:});
            
            t = 'unknownTracer';
            fold = basename(ip.Results.folder);
            if (lstrfind(fold, ip.Results.tracers))
                idx = regexp(fold, ip.Results.prefixRegexp, 'end');
                t = fold(1:idx);
            end
        end
        function v  = visitNumber(str)
            pos = regexp(str, 'V\d');
            assert(~isempty(pos));
            v = str2double(str(pos(end)+1));            
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

