function imgProc_sm(infile, outfile, probename, n, nEvery, projectname, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function is the image processing part of OAP processing using
%  distributed memory parallisation. The function use one simple interface
%  for all probes.
%
%  Interface:
%    infile   :   The input file name
%    outfile  :   The output file name
%    probetype:   One of the following: '2DC','2DP','CIP','PIP','HVPS' and '2DS'
%    n        :   The nth chuck to be processed.
%    nEvery   :   The individual chuck size. nChuck*nEvery shoudl equal the
%                 total frame number
%    projectname: The name of project so that you can write the specific
%                 code for you data
%    varargin :   OPTIONAL argument for the aircraft TAS file name
%                 (HVPS/2DS proves only)
%
%  Note other important variables used in the program
%    handles:  a structure to store information. It is convinient to use a
%          struture to store the global information rather than using
%          various varibles
%
%  Update Dates:
%    * Initially Written by Will Wu, 06/24/2013
%          imgprocdm(File,probetype,n)
%    * Updated by Will Wu, 10/11/2013
%          New function interface
%          imgprocdm(infile,outfile,probetype,n, nEvery) and updated documentation.
%          This version is a major update to include all probes and simplify
%          the function interface significantly
%    * Updated by Will Wu, 07/10/2013
%          New function interface imgProc_dm(infile,outfile,probetype,n, nEvery)
%          Output perimeter, rectangle length/width/angle and eclispe
%          length/width/angle
%    * Added by Wei Wu, May 11, 2016
%          Add the project specific code with projectname in the following format:
%            if strcmp(projectname, 'PECAN')  % For example for PECAN dataset
%               ...
%            end
%    * Updated by Will Wu, 07/11/2016
%          New function name with the option to turn CGAL on and off for
%          speed
%    * Updated by Dan Stechman, 11/28/2016
%          Changed definition of diode_stats variable to sum instances of *shadowed*
%          particles (old version summed illuminated particles).
%          Added boolean and associated code to calculate diode shadow frequency
%          on a particle-by-particle basis.
%    * Updated by Joe Finlon, 03/03/2017
%          Improved calculation of 'Time_in_seconds' variable for HVPS/2DS
%          probes using the TAS from the aircraft data file (still requires
%          diff('Time_in_seconds') to compute int-arr time in
%          post-processing)
%    * Updated by Joe Finlon, 06/22/2017
%          Omit saving DOF/Overload flag for 2DC/2DP since there is no
%          equivalent for these probes.
%          Added metadata for netCDF output.
%    * Updated by Dan Stechman, 08/01/2017
%          Added PECAN-specific PIP time offset correction.
%		   Added [potentially] PECAN-specific CIP/PIP corrupt record identification,
%		   removal, and logging. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setting probe information according to probe type
%    use ProbeType to indicate three type of probes:
%       0: 2DC/2DP, 32 doides, boundary 85,
%       1: CIP/PIP, 64 doides, boundary 170
%       2: HVPS/2DS, 128 doides, boundary 170
	
	
	iRectEllipse = 1;  % Set defualt to no Rectangle fit and Ellipse fit
	
	% Option to save diode stats for every particle
	% More than doubles filesize, and increase computation time
	% Only enable if actually needed
	calcAllDiodeStats = 1;
	
	switch probename
		case '2DC'
			boundary=[255 255 255 255];
			boundarytime=85;
			
			ds = 0.025;			     % Size of diode in millimeters
			handles.diodesize = ds;
			handles.diodenum  = 32;  % Diode number
			handles.current_image = 1;
			probetype=0;
			
		case '2DP'
			boundary=[255 255 255 255];
			boundarytime=85;
			
			ds = 0.200;			     % Size of diode in millimeters
			handles.diodesize = ds;
			handles.diodenum  = 32;  % Diode number
			handles.current_image = 1;
			probetype=0;
			
		case 'CIP'
			boundary=[170, 170, 170, 170, 170, 170, 170, 170];
			boundarytime=NaN;
			
			ds = 0.025;			     % Size of diode in millimeters
			handles.diodesize = ds;
			handles.diodenum  = 64;  % Diode number
			handles.current_image = 1;
			probetype=1;
			
		case 'PIP'
			boundary=[170, 170, 170, 170, 170, 170, 170, 170];
			boundarytime=NaN;
			
			ds = 0.100;			     % Size of diode in millimeters
			handles.diodesize = ds;
			handles.diodenum  = 64;  % Diode number
			handles.current_image = 1;
			probetype=1;
			
		case 'HVPS'
			boundary=[43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690];
			boundarytime=0;
			
			ds = 0.150;			     % Size of diode in millimeters
			handles.diodesize = ds;
			handles.diodenum  = 128; % Diode number
			handles.current_image = 1;
			probetype=2;
			% TAS file info - Added by Joe Finlon - 03/03/17
			tasFile = varargin{1};
			tas = ncread(tasFile, 'TAS'); % aircraft TAS in m/s
			tasTime = ncread(tasFile, 'Time'); % aircraft time in HHMMSS
			disp(['Using TAS data from ', tasFile])
			
		case '2DS'
			boundary=[43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690];
			boundarytime=0;
			
			ds = 0.010;			     % Size of diode in millimeters
			handles.diodesize = ds;
			handles.diodenum  = 128; % Diode number
			handles.current_image = 1;
			probetype=2;
			% TAS file info - Added by Joe Finlon - 03/03/17
			tasFile = varargin{1};
			tas = ncread(tasFile, 'TAS'); % aircraft TAS in m/s
			tasTime = ncread(tasFile, 'Time'); % aircraft time in HHMMSS
			disp(['Using TAS data from ', tasFile])
	end
	
	diodenum = handles.diodenum;
	byteperslice = diodenum/8;
	handles.disagree = 0;
	
	%% Read the particle image files
	handles.f = netcdf.open(infile,'nowrite');
	[~, dimlen] = netcdf.inqDim(handles.f,2);
	[~, handles.img_count] = netcdf.inqDim(handles.f,0);
	size_mat = dimlen;
	warning off all
	diode_stats = zeros(1,diodenum);
	
	
	%% Create output NETCDF file and variables
	f = netcdf.create(outfile, 'clobber');
	dimid0 = netcdf.defDim(f,'time',netcdf.getConstant('NC_UNLIMITED'));
	dimid1 = netcdf.defDim(f,'pos_count',2);
	dimid2 = netcdf.defDim(f,'bin_count',diodenum);
	
	NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
	netcdf.putAtt(f, NC_GLOBAL, 'Software', 'UIOPS/imgProc_sm');
	netcdf.putAtt(f, NC_GLOBAL, 'Institution', 'Univ. Illinois, Dept. Atmos. Sciences');
	netcdf.putAtt(f, NC_GLOBAL, 'Creation Time', datestr(now, 'yyyy/mm/dd HH:MM:SS'));
	netcdf.putAtt(f, NC_GLOBAL, 'Description', ['Contains size, morphological, ',...
		'and other particle properties from an uncompressed image file.']);
	netcdf.putAtt(f, NC_GLOBAL, 'Project', projectname);
	netcdf.putAtt(f, NC_GLOBAL, 'Image Source', infile);
	netcdf.putAtt(f, NC_GLOBAL, 'Probe Type', probename);
	if iRectEllipse && calcAllDiodeStats
		netcdf.putAtt(f, NC_GLOBAL, 'Optional parameters active',...
			'Rectangle & elliptical fits & per-particle diode statistics');
	elseif iRectEllipse && ~calcAllDiodeStats
		netcdf.putAtt(f, NC_GLOBAL, 'Optional parameters active',...
			'Rectangle & elliptical fits');
	elseif ~iRectEllipse && calcAllDiodeStats
		netcdf.putAtt(f, NC_GLOBAL, 'Optional parameters active',...
			'Per-particle diode statistics');
	else
		netcdf.putAtt(f, NC_GLOBAL, 'Optional parameters active', 'None')
	end
	
	
	varid1 = netcdf.defVar(f,'Date','int',dimid0);
	netcdf.putAtt(f, varid1, 'Units', 'YYYYMMDD')
	netcdf.putAtt(f, varid1, 'Description', 'Date of image record in which the particle resides')
	
	varid0  = netcdf.defVar(f,'Time','int',dimid0);
	netcdf.putAtt(f, varid0, 'Units', 'HHMMSS')
	netcdf.putAtt(f, varid0, 'Description', 'Time (UTC) of image record in which the particle resides')
	
	varid2  = netcdf.defVar(f,'msec','short',dimid0);
	netcdf.putAtt(f, varid2, 'Units', 'ms')
	netcdf.putAtt(f, varid2, 'Description', 'Sub-second time of image record in which the particle resides')
	
	varid101  = netcdf.defVar(f,'Time_in_seconds','double',dimid0);
	netcdf.putAtt(f, varid101, 'Units', 'sec')
	netcdf.putAtt(f, varid101, 'Description', 'Time since probe was activated')
	
	if probetype~=0 % Added by Joe Finlon - 06/22/17
		varid102  = netcdf.defVar(f,'SliceCount','short',dimid0);
		netcdf.putAtt(f, varid102, 'Units', '--')
		netcdf.putAtt(f, varid102, 'Description', 'Number of slices containing particle')
		
		varid103  = netcdf.defVar(f,'DMT_DOF_SPEC_OVERLOAD','short',dimid0);
		netcdf.putAtt(f, varid103, 'Units', '--')
		netcdf.putAtt(f, varid103, 'Description', 'Flag denoting out-of-focus (DMT) or overloaded (SPEC) particles')
		
		varid104  = netcdf.defVar(f,'Particle_number_all','int',dimid0);
		netcdf.putAtt(f, varid104, 'Units', '--')
		netcdf.putAtt(f, varid104, 'Description', 'Index of particle in 2-D buffer')
	end
	%varid3 = netcdf.defVar(f,'wkday','double',dimid0);
	varid4  = netcdf.defVar(f,'position','short',[dimid1 dimid0]);
	netcdf.putAtt(f, varid4, 'Units', '--')
	netcdf.putAtt(f, varid4, 'Description', 'Slice within image record where particle is first sampled')
	
	varid5  = netcdf.defVar(f,'particle_time','int',dimid0);
	netcdf.putAtt(f, varid5, 'Units', 'HHMMSS')
	netcdf.putAtt(f, varid5, 'Description', 'Particle time (not available for 2DS/HVPS)')
	
	varid6  = netcdf.defVar(f,'particle_millisec','short',dimid0);
	netcdf.putAtt(f, varid6, 'Units', 'ms')
	netcdf.putAtt(f, varid6, 'Description', 'Sub-second particle time (not available for 2DS/HVPS)')
	
	varid7  = netcdf.defVar(f,'particle_microsec','double',dimid0);
	netcdf.putAtt(f, varid7, 'Units', 'microsec')
	netcdf.putAtt(f, varid7, 'Description', 'Sub-second particle time (unitless clock count for 2DS/HVPS)')
	
	varid8  = netcdf.defVar(f,'parent_rec_num','int',dimid0);
	netcdf.putAtt(f, varid8, 'Units', '--')
	netcdf.putAtt(f, varid8, 'Description', 'Index of image record in which the particle resides')
	
	varid9  = netcdf.defVar(f,'particle_num','short',dimid0);
	netcdf.putAtt(f, varid9, 'Units', '--')
	netcdf.putAtt(f, varid9, 'Description', 'Particle index within current image record')
	
	varid10 = netcdf.defVar(f,'image_length','short',dimid0);
	netcdf.putAtt(f, varid10, 'Units', '--')
	netcdf.putAtt(f, varid10, 'Description', 'Particle length in time direction using # photodiodes')
	
	varid11 = netcdf.defVar(f,'image_width','short',dimid0);
	netcdf.putAtt(f, varid11, 'Units', '--')
	netcdf.putAtt(f, varid11, 'Description', 'Particle length in photodiode direction using # photodiodes')
	
	varid12 = netcdf.defVar(f,'image_area','double',dimid0);
	netcdf.putAtt(f, varid12, 'Units', 'mm^2')
	netcdf.putAtt(f, varid12, 'Description', 'Projected area using the # shadowed photodiodes')
	
	varid13 = netcdf.defVar(f,'image_longest_y','short',dimid0);
	netcdf.putAtt(f, varid13, 'Units', '--')
	netcdf.putAtt(f, varid13, 'Description', 'Longest vertically-oriented chord through particle in time direction')
	
	varid14 = netcdf.defVar(f,'image_max_top_edge_touching','short',dimid0);
	netcdf.putAtt(f, varid14, 'Units', '--')
	netcdf.putAtt(f, varid14, 'Description', 'Maximum # of times the bottom diode is shadowed in succession')
	
	varid15 = netcdf.defVar(f,'image_max_bottom_edge_touching','short',dimid0);
	netcdf.putAtt(f, varid15, 'Units', '--')
	netcdf.putAtt(f, varid15, 'Description', 'Maximum # of times the top diode is shadowed in succession')
	
	varid16 = netcdf.defVar(f,'image_touching_edge','byte',dimid0);
	netcdf.putAtt(f, varid16, 'Units', '--')
	netcdf.putAtt(f, varid16, 'Description', '0 denotes image entirely in array; 1 denotes image touching edge')
	
	varid17 = netcdf.defVar(f,'image_auto_reject','short',dimid0);
	netcdf.putAtt(f, varid17, 'Units', '--')
	netcdf.putAtt(f, varid17, 'Description', ['ASCII reject code ("0" or 48: accepted; "a" or 97: aspect ratio > 6; ',...
		'"t" or 116: aspect ratio > 5 + image touching edge; ',...
		'"p" or 112: < 25% shadowed diodes in rectangle encompassing particle; ',...
		'"h" or 104 or 72 or 117: hollow particle; "s" or 115: split image; ',...
		'"z" or 122: zero area image; "f" or 102: zero area image)'])
	
	varid18 = netcdf.defVar(f,'image_hollow','byte',dimid0);
	netcdf.putAtt(f, varid18, 'Units', '--')
	netcdf.putAtt(f, varid18, 'Description', '0 denotes not hollow; 1 denotes a hollow image')
	
	varid19 = netcdf.defVar(f,'image_center_in','byte',dimid0);
	netcdf.putAtt(f, varid19, 'Units', '--')
	netcdf.putAtt(f, varid19, 'Description', '0 denotes center of particle outside array; 1 denotes center is inside')
	
	varid20 = netcdf.defVar(f,'image_axis_ratio','double',dimid0);
	netcdf.putAtt(f, varid20, 'Units', '--')
	netcdf.putAtt(f, varid20, 'Description', 'Ratio between maximum vertical length and maximum horizontal length')
	
	varid21 = netcdf.defVar(f,'image_diam_circle_fit','double',dimid0);
	netcdf.putAtt(f, varid21, 'Units', 'mm')
	netcdf.putAtt(f, varid21, 'Description', 'Dmax following Heymsfield & Parrish (1978)')
	
	varid22 = netcdf.defVar(f,'image_diam_horiz_chord','double',dimid0);
	netcdf.putAtt(f, varid22, 'Units', 'mm')
	netcdf.putAtt(f, varid22, 'Description', 'Dmax from # slices+1; best for sideways-looking probes')
	
	varid23 = netcdf.defVar(f,'image_diam_horiz_chord_corr','double',dimid0);
	netcdf.putAtt(f, varid23, 'Units', 'mm')
	netcdf.putAtt(f, varid23, 'Description', 'Dmax from # slices+1, with Korolev (2007) correction applied to hollow spherical particles')
	
	varid24 = netcdf.defVar(f,'image_diam_following_bamex_code','double',dimid0);
	netcdf.putAtt(f, varid24, 'Units', 'mm')
	netcdf.putAtt(f, varid24, 'Description', 'Dmax from maximum length chord through the particle')
	
	varid25 = netcdf.defVar(f,'image_diam_vert_chord','double',dimid0);
	netcdf.putAtt(f, varid25, 'Units', 'mm')
	netcdf.putAtt(f, varid25, 'Description', ['Dmax from maximum length in photodiode direction; ',...
		'best for sideways-looking probes'])
	
	varid26 = netcdf.defVar(f,'image_diam_minR','double',dimid0);
	netcdf.putAtt(f, varid26, 'Units', 'mm')
	netcdf.putAtt(f, varid26, 'Description', 'Dmax of smallest circle enclosing the particle')
	
	varid27 = netcdf.defVar(f,'image_diam_AreaR','double',dimid0);
	netcdf.putAtt(f, varid27, 'Units', 'mm')
	netcdf.putAtt(f, varid27, 'Description', 'Area equivalent diameter')
	
	varid45 = netcdf.defVar(f,'image_perimeter','double',dimid0);
	netcdf.putAtt(f, varid45, 'Units', 'mm')
	netcdf.putAtt(f, varid45, 'Description', 'Perimeter following the particle boundary')
	
	if 1==iRectEllipse
		varid46 = netcdf.defVar(f,'image_RectangleL','double',dimid0);
		netcdf.putAtt(f, varid46, 'Units', 'mm')
		netcdf.putAtt(f, varid46, 'Description', 'Length of smallest rectangle enclosing the particle')
		
		varid47 = netcdf.defVar(f,'image_RectangleW','double',dimid0);
		netcdf.putAtt(f, varid47, 'Units', 'mm')
		netcdf.putAtt(f, varid47, 'Description', 'Width of smallest rectangle enclosing the particle')
		
		varid67 = netcdf.defVar(f,'image_RectangleAngle','double',dimid0);
		netcdf.putAtt(f, varid67, 'Units', 'radians')
		netcdf.putAtt(f, varid67, 'Description', 'Angle of rectangle with respect to the time direction')
		
		varid48 = netcdf.defVar(f,'image_EllipseL','double',dimid0);
		netcdf.putAtt(f, varid48, 'Units', 'mm')
		netcdf.putAtt(f, varid48, 'Description', 'Length of smallest ellipse enclosing the particle')
		
		varid49 = netcdf.defVar(f,'image_EllipseW','double',dimid0);
		netcdf.putAtt(f, varid49, 'Units', 'mm')
		netcdf.putAtt(f, varid49, 'Description', 'Width of smallest ellipse enclosing the particle')
		
		varid69 = netcdf.defVar(f,'image_EllipseAngle','double',dimid0);
		netcdf.putAtt(f, varid69, 'Units', 'radians')
		netcdf.putAtt(f, varid69, 'Description', 'Angle of rectangle with respect to the time direction')
	end
	varid28 = netcdf.defVar(f,'percent_shadow_area','double',dimid0);
	netcdf.putAtt(f, varid28, 'Units', 'percent')
	netcdf.putAtt(f, varid28, 'Description', 'Ratio between the projected area and L*W')
	
	varid29 = netcdf.defVar(f,'edge_at_max_hole','short',dimid0);
	netcdf.putAtt(f, varid29, 'Units', '--')
	netcdf.putAtt(f, varid29, 'Description', ['# diodes between edges of the particle for the slice ',...
		'containing the largest gap inside the particle'])
	
	varid30 = netcdf.defVar(f,'max_hole_diameter','short',dimid0);
	netcdf.putAtt(f, varid30, 'Units', '--')
	netcdf.putAtt(f, varid30, 'Description', 'Diameter of the largest hole inside the particle')
	
	varid31 = netcdf.defVar(f,'part_z','double',dimid0);
	netcdf.putAtt(f, varid31, 'Units', '--')
	netcdf.putAtt(f, varid31, 'Description', 'Particle depth in object plane calculated via lookup table')
	
	varid32 = netcdf.defVar(f,'size_factor','double',dimid0);
	netcdf.putAtt(f, varid32, 'Units', '--')
	netcdf.putAtt(f, varid32, 'Description', 'Dmax reduction factor folowing Korolev (2007) correction')
	
	varid33 = netcdf.defVar(f,'holroyd_habit','short',dimid0);
	netcdf.putAtt(f, varid33, 'Units', '--')
	netcdf.putAtt(f, varid33, 'Description', ['ASCII habit code following the Holroyd (1987) algorithm ',...
		'("M" or 77: zero image; "C" or 67: center-out image; "t" or 116: tiny; "o" or 111: oriented; ',...
		'"l" or 108: linear; "a" or 97: aggregate; "g" or 103: graupel; "s" or 115: sphere; ',...
		'"h" or 104: hexagonal; "i" or 105: irregular; "d" or 100: dendrite)'])
	
	varid34 = netcdf.defVar(f,'area_hole_ratio','double',dimid0);
	netcdf.putAtt(f, varid34, 'Units', '--')
	netcdf.putAtt(f, varid34, 'Description', 'Ratio between the projected area and the area of hole inside particle')
	
	varid35 = netcdf.defVar(f,'inter_arrival','double',dimid0);
	netcdf.putAtt(f, varid35, 'Units', 'sec')
	netcdf.putAtt(f, varid35, 'Description', ['Inter-arrival time between particles ',...
		'(use diff(time_in_seconds) for 2DS/HVPS)'])
	
	varid36 = netcdf.defVar(f,'bin_stats','int',dimid2);
	netcdf.putAtt(f, varid36, 'Units', '--')
	netcdf.putAtt(f, varid36, 'Description', '# times specified photodiode is shadowed for particles in this file')
	
	if calcAllDiodeStats
		varid37 = netcdf.defVar(f,'image_bin_stats','int',[dimid2 dimid0]);
		netcdf.putAtt(f, varid37, 'Units', '--')
		netcdf.putAtt(f, varid37, 'Description', '# times specified photodiode is shadowed for each particle')
	end
	
	netcdf.endDef(f)
	
	%% Variables initialization
	kk=1;
	w=-1;
	wstart = 0;
	
	time_offset_hr = 0;
	time_offset_mn = 0;
	time_offset_sec = 0;
	time_offset_ms = 0;
	timeset_flag = 0;
	
	% Initialize corrupt record and particle counters - Added by Dan Stechman 8/1/17
	% Currently, issue is specific to PECAN, but may be present with other DMT data
	if probetype == 1
		crptRecCount = 0;
		crptPartCount = 0;
	end
	
	%% Processing nth chuck. Every chuck is nEvery frame
	%% Analyze each individual particle images and Outpu the particle by particle information
	for i=((n-1)*nEvery+1):min(n*nEvery,handles.img_count)
		
		handles.year     = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'year'    ),i-1,1);
		handles.month    = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'month'  ),i-1,1);
		handles.day      = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'day'  ),i-1,1);
		handles.hour     = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'hour'    ),i-1,1);
		handles.minute   = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'minute'  ),i-1,1);
		handles.second   = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'second'  ),i-1,1);
		handles.millisec = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'millisec'),i-1,1);
		
		if mod(i,100) == 0
			disp([num2str(i),'/',num2str(handles.img_count), ', ',datestr(now)])
		end
		varid = netcdf.inqVarID(handles.f,'data');
		
		if probetype==0
			temp = netcdf.getVar(handles.f,varid,[0, 0, i-1], [4,1024,1]);
		else
			temp = netcdf.getVar(handles.f,varid,[0, 0, i-1], [8,1700,1]);
		end
		data(:,:) = temp';
		
		j=1;
		start=0;
		firstpart = 1;
		
		% Check the current data record for corrupt boundaries and keep track
		% of which and how many records and particles are being dicarded
		% May be PECAN-specific, but should be tested with other DMT probe data - Added by Dan Stechman 8/1/17
		if probetype == 1
			tempMmbr = ismember(data,boundary);
			memberSum = sum(tempMmbr,2);
			numCrptBnds = length(find(memberSum > 4 & memberSum < 8));
			numGoodBnds = length(find(memberSum == 8));

			if numCrptBnds > 0
				fprintf('%d/%d boundaries in record %d are corrupt.\n',numCrptBnds,numGoodBnds+numCrptBnds,i);
				crptRecCount = crptRecCount + 1;
				crptPartCount = crptPartCount + (numGoodBnds + numCrptBnds);
			end
		end
		
		%c=[dec2bin(data(:,1),8),dec2bin(data(:,2),8),dec2bin(data(:,3),8),dec2bin(data(:,4),8)];
		while data(j,1) ~= -1 && j < size(data,1) && numCrptBnds < 1
			% Calculate every particles
			if (isequal(data(j,:), boundary) && ( (isequal(data(j+1,1), boundarytime) || probetype==1) ) )
				if start == 0
					if probetype == 0 || probetype == 1
						start = 2;
					else
						start = 1;
					end
				end
				
				if probetype == 0
					if start+1 > (j-1)  % Remove Corrupted Data
						break;
					end
				else
					if start > (j-1)  % Remove Corrupted Data
						break;
					end
				end
				
				header_loc = j+1;
				w=w+1;
				%% Create binary image according to probe type
				
				if probetype==0
					ind_matrix(1:j-start-1,:) = data(start+1:j-1,:);  % 2DC has 3 slices between particles (sync word timing word and end of particle words)
					c=[dec2bin(ind_matrix(:,1),8),dec2bin(ind_matrix(:,2),8),dec2bin(ind_matrix(:,3),8),dec2bin(ind_matrix(:,4),8)];
				elseif probetype==1
					ind_matrix(1:j-start,:) = data(start:j-1,:);
					c=[dec2bin(ind_matrix(:,1),8), dec2bin(ind_matrix(:,2),8),dec2bin(ind_matrix(:,3),8),dec2bin(ind_matrix(:,4),8), ...
						dec2bin(ind_matrix(:,5),8), dec2bin(ind_matrix(:,6),8),dec2bin(ind_matrix(:,7),8),dec2bin(ind_matrix(:,8),8)];
				elseif probetype==2
					ind_matrix(1:j-start,:) = 65535 - data(start:j-1,:); % I used 1 to indicate the illuminated doides for HVPS
					c=[dec2bin(ind_matrix(:,1),16), dec2bin(ind_matrix(:,2),16),dec2bin(ind_matrix(:,3),16),dec2bin(ind_matrix(:,4),16), ...
						dec2bin(ind_matrix(:,5),16), dec2bin(ind_matrix(:,6),16),dec2bin(ind_matrix(:,7),16),dec2bin(ind_matrix(:,8),16)];
				end
				
				% Just to test if there are bad images, usually 0 area images
				figsize = size(c);
				if figsize(2)~=diodenum
					disp('Not equal to doide number');
					return
				end
				
				
				images.position(kk,:) = [start, j-1];
				parent_rec_num(kk)=i;
				particle_num(kk) = mod(kk,66536); %hex2dec([dec2hex(data(start-1,7)),dec2hex(data(start-1,8))]);
				
				%  Get the particle time
				if probetype == 0
					bin_convert = [dec2bin(data(header_loc,2),8),dec2bin(data(header_loc,3),8),dec2bin(data(header_loc,4),8)];
					part_time = bin2dec(bin_convert);       % Interarrival time in tas clock cycles
					tas2d = netcdf.getVar(handles.f,netcdf.inqVarID(handles.f,'tas'),i-1, 1);
					part_time = part_time/tas2d*handles.diodesize/(10^3);
					time_in_seconds(kk) = part_time;
					
					images.int_arrival(kk) = part_time;
					
					if(firstpart == 1)
						firstpart = 0;
						start_hour = handles.hour;
						start_minute = handles.minute;
						start_second = handles.second;
						start_msec = handles.millisec*10;
						% First, we get the hours....
						start_msec = start_msec;
						start_microsec = 0;
						time_offset_hr = 0;
						time_offset_mn = 0;
						time_offset_sec = 0;
						time_offset_ms = 0;
						
						part_hour(kk) = start_hour;
						part_min(kk) = start_minute;
						part_sec(kk) = start_second;
						part_mil(kk) = start_msec;
						part_micro(kk) = 0;
					else
						frac_time = part_time - floor(part_time);
						frac_time = frac_time * 1000;
						part_micro(kk) = part_micro(kk-1) + (frac_time - floor(frac_time))*1000;
						part_mil(kk) = part_mil(kk-1) + floor(frac_time);
						part_sec(kk) = part_sec(kk-1) + floor(part_time);
						part_min(kk) = part_min(kk-1);
						part_hour(kk) = part_hour(kk-1);
					end
					
					part_sec(part_mil >= 1000) = part_sec(part_mil >= 1000) + 1;
					part_mil(part_mil >= 1000) = part_mil(part_mil >= 1000) - 1000;
					
					part_min(part_sec >= 60) = part_min(part_sec >= 60) + 1;
					part_sec(part_sec >= 60) = part_sec(part_sec >= 60) - 60;
					
					part_hour(part_min >= 60) = part_hour(part_min >= 60) + 1;
					part_min(part_min >= 60) = part_min(part_min >= 60) - 60;
					part_hour(part_hour >= 24) = part_hour(part_hour >= 24) - 24;
				
				elseif probetype == 1
					bin_convert = [dec2bin(data(start-1,2),8),dec2bin(data(start-1,3),8),dec2bin(data(start-1,4),8), ...
						dec2bin(data(start-1,5),8), dec2bin(data(start-1,6),8)];
					
					part_hour(kk) = bin2dec(bin_convert(1:5));
					part_min(kk) = bin2dec(bin_convert(6:11));
					part_sec(kk) = bin2dec(bin_convert(12:17));
					part_mil(kk) = bin2dec(bin_convert(18:27));
					part_micro(kk) = bin2dec(bin_convert(28:40))*125e-9;
					
					% Apply time correction to PIP probe-time for slect PECAN flights - Added by Dan Stechman 8/1/17
					if strcmp(projectname, 'PECAN') && strcmp(probename, 'PIP')
						if (handles.month == 6 && handles.day == 17)
							part_hour_offset = -12;
							part_min_offset = 0;
							part_sec_offset = 0;
						elseif (handles.month == 6 && handles.day == 20) % Time offset was corrected during flight after 6 UTC
							if handles.hour < 6
								part_hour_offset = -12;
								part_min_offset = 0;
								part_sec_offset = 0;
							else
								part_hour_offset = 0;
								part_min_offset = 0;
								part_sec_offset = 0;
							end
						elseif (handles.month == 7 && (handles.day >= 6 && handles.day <= 13))
							part_hour_offset = -12;
							part_min_offset = 0;
							part_sec_offset = 0;
						else
							part_hour_offset = 0;
							part_min_offset = 0;
							part_sec_offset = 0;
						end
						
						
						part_hour(kk) = part_hour(kk) + part_hour_offset;
						part_min(kk) = part_min(kk) + part_min_offset;
						part_sec(kk) = part_sec(kk) + part_sec_offset;
						
					end
					
					particle_sliceCount(kk)=bitand(data(start-1,1),127);
					particle_DOF(kk)=bitand(data(start-1,1),128);
					particle_partNum(kk)=bin2dec([dec2bin(data(start-1,7),8),dec2bin(data(start-1,8),8)]);
					
					time_in_seconds(kk) = part_hour(kk) * 3600 + part_min(kk) * 60 + part_sec(kk) + part_mil(kk)/1000 + part_micro(kk);
					
					
					if kk > 1
						images.int_arrival(kk) = time_in_seconds(kk) - time_in_seconds(kk-1);
					else
						images.int_arrival(kk) = 0;
					end
					
					
				elseif probetype==2
					
					particle_DOF(kk)=bitand(data(header_loc,4), 32768);
					particle_partNum(kk)=double(data(header_loc,5));
					particle_sliceCount(kk)=double(data(header_loc,6));
					
					part_time = double(data(header_loc,7))*2^16+double(data(header_loc,8));       % Interarrival time in clock cycles
					part_micro(kk) = part_time;
					part_mil(kk)   = 0;
					part_sec(kk)   = 0;
					part_min(kk)   = 0;
					part_hour(kk)  = 0;
					
					% -----------------------------------------------------
					% calculate the particle time (in seconds) - Added by Joe Finlon - 03/03/17
					tempTime = double(handles.hour)*10000+double(handles.minute)*100+double(handles.second); % image record time in HHMMSS
					tas2d = tas(tasTime==tempTime); % get TAS for corresponding time period
					
					if exist('part_time_prev', 'var')
						if isempty(tas2d) % use when no TAS data is available for current time
							% compute time difference (sec) between
							% particles to determine whether corrections
							% are needed (see 'indexRollback' in
							% sizeDist.m)
							timeDiff = (part_time-part_time_prev)*(handles.diodesize/(10^3)/170);
							
							if timeDiff>=-250 % clock cycle hasn't rolled back
								time_in_seconds(kk) = time_in_seconds_prev + timeDiff;
							else % perform correction when clock cycle exceeds 2^32 and rolls back
								time_in_seconds(kk) = time_in_seconds_prev + (2^32-1-part_time_prev+part_time)*...
									(handles.diodesize/(10^3)/170);
								disp(['Performing clock cyle time correction. ',num2str(part_time_prev),...
									' ',num2str(part_time),' ',num2str(time_in_seconds_prev),' ',num2str(time_in_seconds(kk))])
							end
						else % use the TAS to convert clock cycles to time
							% compute time difference (sec) between
							% particles to determine whether corrections
							% are needed (see 'indexRollback' in
							% sizeDist.m)
							timeDiff = (part_time-part_time_prev)*(handles.diodesize/(10^3)/tas2d);
							
							if timeDiff>=-250 % clock cycle hasn't rolled back
								time_in_seconds(kk) = time_in_seconds_prev + timeDiff;
							else % perform correction when TAS clock cycle exceeds 2^32 and rolls back
								time_in_seconds(kk) = time_in_seconds_prev + (2^32-1-part_time_prev+part_time)*...
									(handles.diodesize/(10^3)/tas2d);
								disp(['Performing clock cyle time correction. ',num2str(part_time_prev),...
									' ',num2str(part_time),' ',num2str(time_in_seconds_prev),' ',num2str(time_in_seconds(kk))])
							end
						end
					else % run for the first particle in data file
						if isempty(tas2d)
							time_in_seconds(kk) = part_time*(handles.diodesize/(10^3)/170);
						else
							time_in_seconds(kk) = part_time*(handles.diodesize/(10^3)/tas2d);
						end
					end
					
					if(kk>1) % ** Use diff(time_in_seconds) to compute int-arr time in post-processing **
						images.int_arrival(kk) = time_in_seconds(kk)-time_in_seconds(kk-1);
					else
						images.int_arrival(kk) = 0;
					end
					
					part_time_prev = part_time; % assign the clock cycle for use in next iteration
					time_in_seconds_prev = time_in_seconds(kk);
					% -----------------------------------------------------
				end
				
				temptimeinhhmmss = part_hour(kk) * 10000 + part_min(kk) * 100 + part_sec(kk);
				%                 if (temptimeinhhmmss<0 || temptimeinhhmmss>240000)
				%                    fprintf('%d\n',temptimeinhhmmss);
				%                 end
				
				slices_ver = length(start:j-1);
				rec_time(kk)=double(handles.hour)*10000+double(handles.minute)*100+double(handles.second);
				rec_date(kk)=double(handles.year)*10000+double(handles.month)*100+double(handles.day);
				rec_millisec(kk)=handles.millisec;
				%                 rec_wkday(kk)=handles.wkday(i);
				
				
				% % 				if mod(i,100) == 0
				% % 					fprintf('\trecord time - particle time = %.3f\n',(hhmmss2insec(rec_time(kk))-time_in_seconds(kk)));
				% % 				end
				
				
				%% Determine the Particle Habit
				%  We use the Holroyd algorithm here
				handles.bits_per_slice = diodenum;
				if calcAllDiodeStats
					images.diode_stats(kk,:) = sum(c=='0',1); % Option to save diode stats for every particle
				end
				diode_stats = diode_stats + sum(c=='0',1);
				csum = sum(c=='0',1);
				
				images.holroyd_habit(kk) = holroyd(handles,c);
				
				%% Determine if the particle is rejected or not
				%  Calculate the Particle Length, Width, Area, Auto Reject
				%  Status And more... See calculate_reject_unified()
				%  funtion for more information
				
				[images.image_length(kk),images.image_width(kk),images.image_area(kk), ...
					images.longest_y_within_a_slice(kk),images.max_top_edge_touching(kk),images.max_bottom_edge_touching(kk),...
					images.image_touching_edge(kk), images.auto_reject(kk),images.is_hollow(kk),images.percent_shadow(kk),images.part_z(kk),...
					images.sf(kk),images.area_hole_ratio(kk),handles]=calculate_reject_unified(c,handles,images.holroyd_habit(kk));
				
				images.max_hole_diameter(kk) = handles.max_hole_diameter;
				images.edge_at_max_hole(kk) = handles.edge_at_max_hole;
				
				max_horizontal_length = images.image_length(kk);
				max_vertical_length = images.longest_y_within_a_slice(kk);
				image_area = images.image_area(kk);
				
				diode_size= handles.diodesize;
				corrected_horizontal_diode_size = handles.diodesize;
				largest_edge_touching  = max(images.max_top_edge_touching(kk), images.max_bottom_edge_touching(kk));
				smallest_edge_touching = min(images.max_top_edge_touching(kk), images.max_bottom_edge_touching(kk));
				
				%% Calculate more size deciptor using more advanced techniques
				%  See dropsize for more information
				[images.center_in(kk),images.axis_ratio(kk),images.diam_circle_fit(kk),images.diam_horiz_chord(kk),images.diam_vert_chord(kk),...
					images.diam_horiz_mean(kk), images.diam_spheroid(kk)]=dropsize(max_horizontal_length,max_vertical_length,image_area...
					,largest_edge_touching,smallest_edge_touching,diode_size,corrected_horizontal_diode_size, diodenum);
				
				%% Calculate size deciptor using bamex code
				%  See dropsize_new for more information
				% images.diam_bamex(kk) = dropsize_new(c, largest_edge_touching, smallest_edge_touching, diodenum, corrected_horizontal_diode_size, handles.diodesize, max_vertical_length);
				
				%% Using OpenCV C program to calculate length, width and radius. This
				%% Get diameter of the smallest-enclosing circle, rectangle and ellipse
				%images.minR(kk)=particlesize_cgal(c);
				images.minR(kk)=CGAL_minR(c);
				images.AreaR(kk)=2*sqrt(images.image_area(kk)/3.1415926);  % Calculate the Darea (area-equivalent diameter)
				images.Perimeter(kk)=ParticlePerimeter(c);
				
				if 1==iRectEllipse
					[images.RectangleL(kk), images.RectangleW(kk), images.RectangleAngle(kk)] = CGAL_RectSize(c);
					[images.EllipseL(kk), images.EllipseW(kk), images.EllipseAngle(kk)]       = CGAL_EllipseSize(c);
				end
				%% Get the area ratio using the DL=max(DT,DP), only observed area are used
				if images.image_length(kk) > images.image_width(kk)
					images.percent_shadow(kk) = images.image_area(kk) / (pi * images.image_length(kk).^ 2 / 4);
				elseif images.image_width(kk) ~= 0
					images.percent_shadow(kk) = images.image_area(kk) / (pi * images.image_width(kk).^ 2 / 4);
				else
					images.percent_shadow(kk) = 0;
				end
				
				start = j + 2;
				kk = kk + 1;
				clear c ind_matrix
			end
			
			j = j + 1;
		end
		
		%% Write out the processed information on NETCDF
		if kk > 1
			
			
			netcdf.putVar ( f, varid0, wstart, w-wstart+1, rec_time(:) );
			netcdf.putVar ( f, varid1, wstart, w-wstart+1, rec_date(:) );
			
			netcdf.putVar ( f, varid101, wstart, w-wstart+1, time_in_seconds(:) );
			if probetype~=0 % Added by Joe Finlon - 06/22/17
				netcdf.putVar ( f, varid102, wstart, w-wstart+1, particle_sliceCount );
				netcdf.putVar ( f, varid103, wstart, w-wstart+1, particle_DOF );
				netcdf.putVar ( f, varid104, wstart, w-wstart+1, particle_partNum );
			end
			
			netcdf.putVar ( f, varid2, wstart, w-wstart+1, rec_millisec(:) );
			%netcdf.putVar ( f, varid3, wstart, w-wstart+1, rec_wkday(:) );
			netcdf.putVar ( f, varid4, [0 wstart], [2 w-wstart+1], images.position' );
			netcdf.putVar ( f, varid5, wstart, w-wstart+1, part_hour(:)*10000+part_min(:)*100+part_sec(:) );
			netcdf.putVar ( f, varid6, wstart, w-wstart+1, part_mil(:) );
			netcdf.putVar ( f, varid7, wstart, w-wstart+1, part_micro(:) );
			netcdf.putVar ( f, varid8, wstart, w-wstart+1, parent_rec_num );
			netcdf.putVar ( f, varid9, wstart, w-wstart+1, particle_num(:) );
			netcdf.putVar ( f, varid10, wstart, w-wstart+1, images.image_length);
			netcdf.putVar ( f, varid11, wstart, w-wstart+1, images.image_width);
			netcdf.putVar ( f, varid12, wstart, w-wstart+1, images.image_area*diode_size*diode_size);
			netcdf.putVar ( f, varid13, wstart, w-wstart+1, images.longest_y_within_a_slice);
			netcdf.putVar ( f, varid14, wstart, w-wstart+1, images.max_top_edge_touching);
			netcdf.putVar ( f, varid15, wstart, w-wstart+1, images.max_bottom_edge_touching);
			netcdf.putVar ( f, varid16, wstart, w-wstart+1, images.image_touching_edge-'0');
			netcdf.putVar ( f, varid17, wstart, w-wstart+1, int16(images.auto_reject));
			netcdf.putVar ( f, varid18, wstart, w-wstart+1, images.is_hollow);
			netcdf.putVar ( f, varid19, wstart, w-wstart+1, images.center_in);
			netcdf.putVar ( f, varid20, wstart, w-wstart+1, images.axis_ratio);
			netcdf.putVar ( f, varid21, wstart, w-wstart+1, images.diam_circle_fit);
			netcdf.putVar ( f, varid22, wstart, w-wstart+1, images.diam_horiz_chord);
			netcdf.putVar ( f, varid23, wstart, w-wstart+1, images.diam_horiz_chord ./ images.sf);
			netcdf.putVar ( f, varid24, wstart, w-wstart+1, images.diam_horiz_mean);
			netcdf.putVar ( f, varid25, wstart, w-wstart+1, images.diam_vert_chord);
			netcdf.putVar ( f, varid26, wstart, w-wstart+1, images.minR*diode_size);
			netcdf.putVar ( f, varid27, wstart, w-wstart+1, images.AreaR*diode_size);
			netcdf.putVar ( f, varid45, wstart, w-wstart+1, images.Perimeter*diode_size);
			if 1==iRectEllipse
				netcdf.putVar ( f, varid46, wstart, w-wstart+1, images.RectangleL*diode_size);
				netcdf.putVar ( f, varid47, wstart, w-wstart+1, images.RectangleW*diode_size);
				netcdf.putVar ( f, varid67, wstart, w-wstart+1, images.RectangleAngle);
				netcdf.putVar ( f, varid48, wstart, w-wstart+1, images.EllipseL*diode_size);
				netcdf.putVar ( f, varid49, wstart, w-wstart+1, images.EllipseW*diode_size);
				netcdf.putVar ( f, varid69, wstart, w-wstart+1, images.EllipseAngle);
			end
			netcdf.putVar ( f, varid28, wstart, w-wstart+1, images.percent_shadow);
			netcdf.putVar ( f, varid29, wstart, w-wstart+1, images.max_hole_diameter);
			netcdf.putVar ( f, varid30, wstart, w-wstart+1, images.edge_at_max_hole);
			netcdf.putVar ( f, varid31, wstart, w-wstart+1, images.part_z);
			netcdf.putVar ( f, varid32, wstart, w-wstart+1, images.sf);
			netcdf.putVar ( f, varid33, wstart, w-wstart+1, int16(images.holroyd_habit));
			netcdf.putVar ( f, varid34, wstart, w-wstart+1, images.area_hole_ratio);
			netcdf.putVar ( f, varid35, wstart, w-wstart+1, images.int_arrival);
			netcdf.putVar ( f, varid36, diode_stats );
			
			if calcAllDiodeStats
				netcdf.putVar ( f, varid37, [0 wstart], [diodenum w-wstart+1], images.diode_stats' ); % Option to save diode stats for every particle
			end
			
			wstart = w+1;
			kk = 1;
			clear rec_time rec_date rec_millisec part_hour part_min part_sec part_mil part_micro parent_rec_num particle_num images time_in_seconds particle_sliceCount particle_DOF particle_partNum
			
		end
		clear images
		
	end
	
	if probetype == 1
		slicesInChunk = length(((n-1)*nEvery+1):min(n*nEvery,handles.img_count));
		fprintf('Chunk %d complete -- %d/%d records were likely corrupt (%.4f%%) --> %d potential particles ignored.\n',n,crptRecCount,...
			slicesInChunk,(crptRecCount/slicesInChunk),crptPartCount);
	end
	
	warning on all
	
	netcdf.close(f);
end