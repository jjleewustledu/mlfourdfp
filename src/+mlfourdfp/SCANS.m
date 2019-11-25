classdef SCANS 
	%% SCANS  

	%  $Revision$
 	%  was created 29-Oct-2019 17:06:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlfourdfp/src/+mlfourdfp.
 	%% It was developed on Matlab 9.7.0.1216025 (R2019b) Update 1 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	properties
        mrdir
    end
    
    properties (Dependent)
        description
        lengthSeries
        sequence
        series
    end

	methods (Static)
        function this = buildScansStudiesTxt(mrdir)
            pwd0 = pushd(mrdir);
            system(fullfile(getenv('RELEASE'), 'pseudo_dcm_sort.csh SCANS -s'))
            this = mlfourdfp.SCANS(mrdir);
            popd(pwd0)
        end
    end
    
    methods
        
        %% GET
        
        function g = get.description(this)
            g = this.asTable_.description;
        end
        function g = get.lengthSeries(this)
            g = this.asTable_.lengthSeries;
        end
        function g = get.sequence(this)
            g = this.asTable_.sequence;
        end
        function g = get.series(this)
            g = this.asTable_.series;
        end
        
        
        %%
        
        function tbl = table(this)
            tbl = this.asTable_;
        end
		  
 		function this = SCANS(varargin)
 			%% SCANS builds SCANS.studies.txt as needed and stores its contents as a table.
 			%  @param mrdir.

            ip = inputParser;
            addRequied(ip, 'mrdir', @isfolder)
            parse(ip, varargin{:})
            this.mrdir = ip.Results.mrdir;

            if ~isfile(fullfile(this.mrdir, 'SCANS.studies.txt'))
                this.buildScansStudiesTxt(this.mrdir);
            end
            opts = detectImportOptions(fullfile(this.mrdir, 'SCANS.studies.txt'));
            opts.VariableNames = {'series' 'sequence' 'description' 'lengthSeries'};
            this.asTable_ = readtable(fullfile(this.mrdir, 'SCANS.studies.txt'), opts);
        end
    end
    
    %% PRIVATE
    
    properties (Access = private)
        asTable_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

