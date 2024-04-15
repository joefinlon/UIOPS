% holroyd - identified particle habit according to Holroyd  (1987)
% inputs:
%    handles - handles structure outlined in run_img_processing.m
%    image_buffer - n x photodiodes/8 raw image buffer without timestamps
% outputs:
%    holroyd_habit - habit code as listed below
function [holroyd_habit, F, w_l_p, c_r, e_d_r, m_g, O, a_r, p] = holroyd_new(handles, image_buffer)

%/***************************************************************/%
%/*	Return code                                                */
%/*	                                                           */
%/*	                                                           */
%/*	references: 											   */
%/* original methodology: J. Atmos. and Oceanic Tech. Vol 4,   */
%/* Sept. '87 pages 498- 511. 	                               */
%/*	modified according as outlined in Schima et al. 2024       */
%/*    (A Multi-Probe Automated Classification of Ice Crystal  */
%/*    Habits During the IMPACTS Campaign)                     */
%/*	                                                           */
%/*	'M' = not calculated, zero image                           */
%/*	'C' = not calculated, center is out                        */
%/*	't' = tiny                                                 */
%/*	'o' = bad                                                  */
%/*	'l' = linear                                               */
%/*	'a' = aggregate                                            */
%/*	'g' = graupel                                              */
%/*	's' = spherical                                            */
%/*	'h' = hexagonal                                            */
%/*	'i' = irregular                                            */
%/*	'd' = dendrite                                             */
%/*	                                                           */
%/***************************************************************/

	image_size = size(image_buffer);
    probe_resolution = handles.diodesize;
	n_slices  = image_size(1);
	
	% Set variables to -1 initially; this indicates that they are not valid if they are still below zero after processing
	F = -1;
	w_l_p = -1;
	c_r = -1;
	e_d_r = -1;
	m_g = -1;
	O = -1;
	a_r = -1;
		
	% Tuning coefficients; used to adjust algorithm for resolution differences
	
	% Sphere coefficent; higher = easier to get spheres
	C_1 = 1.0;
	
	% Graupel coefficient; higher = easier to get graupel
	C_2 = 1.0;
	
	% Dendrite coefficient; higher = harder to get dendrites
	C_3 = 1.0;
	
	% Aggregate coefficient; higher = harder to get aggregates
	C_4 = 1.0;
	
	% Plate coefficient; higher = harder to get plates
	C_5 = 1.0;

	if (n_slices == 0)
	    holroyd_habit = 'M';
	    return;
	else
		if (parabola_fit_center_is_in(image_buffer, n_slices) == 1) 
			[x, y, d, a_a, a, r_2, F, O, w_l_p, w_l_m_p, s_c, c_r, e_d_r, m_g, a_r, p] = calc_stat(handles,image_buffer, n_slices);
			
			% Convert measurements to microns; assumes probe resolution given in mm
			d = d * probe_resolution * 1000;
			x = x * probe_resolution * 1000;
			y = y * probe_resolution * 1000;
			%d = d * probe_resolution * 1000;
			a = a * (probe_resolution * 1000) * (probe_resolution * 1000);
			w_l_p =  w_l_p * probe_resolution * 1000;
			w_l_m_p = w_l_m_p * probe_resolution * 1000;
			m_g = m_g * probe_resolution * 1000;
			p = p * probe_resolution * 1000;

			% Missing image
			if (a == 0 )
				holroyd_habit = 'M';
				return;
			end
	
			% Tiny particle (too small to classify; what is too small depends on probe resolution)
			if (d < 300)
				holroyd_habit = 't';
				return;
			end
			
			% Removal of bad particles with large internal gaps 
			if ((m_g > (d/3)) || (m_g >= 10 * (probe_resolution * 1000))) || (a_r > 20)
				holroyd_habit = 'o';
				return;
			end
			
			if e_d_r > 0.5   
            	% Assume particle is center-out because too much is touching the edge
        		holroyd_habit = 'C';
        		return;
        	end
			
			% Calculate required aspect ratio to define a particle as a column
			needed_ratio = (1 + ((w_l_p * w_l_m_p)/40000) + (60/(d * (probe_resolution * 1000)))*(w_l_p) + (d * 0.0005));
			
			% An aspect ratio of 5 is always sufficient for a column
			if (needed_ratio > 5)
            	needed_ratio = 5;
            end
			
			if (a_r >= needed_ratio)
				% Column
				holroyd_habit = 'l';
				return;

			end
			
			% Calculate needed fine-detail ratios for certain habit classes; this algorithm
			% heavily relies on the fine detail ratio for calculations.
			sphere_f = C_1 * ((13 - (d / 250) - (abs(1 - c_r)*7.5)) / a_r);
			
			if d < 1000
				graupel_f = C_2 * (25 - ((1 - O) * 15) + (d - 1000)/30 - (abs(1 - c_r)*15)) / a_r;
			else
				graupel_f = C_2 * (25 - ((1 - O) * 15) + (d - 1000)/300 - (abs(1 - c_r)*15)) / a_r;
			end

			if d < 1000
				dendrite_f = C_3 * (35 + d*d/3e5) * a_r * (1 + (1000 - d)/200);
			else
				dendrite_f = C_3 * (35 + d*d/3e5) * a_r;
			end
			
			if d < 1000
				aggregate_f = C_4 * (30 - (d - 1000) / 30 + s_c * 10);
			else
				aggregate_f = C_4 * (30 - (d - 1000) / 300 + s_c * 10);
			end
			
			plate_f = C_5 * (50 - (((60 - w_l_m_p) / 10) * 8 / a_r)) * (2 * e_d_r + 1) * (1 + (2500/(d * d)) * probe_resolution);
			if d > 1000
				plate_f = plate_f * (d / 1000);
			end

			% Habit is a sphere
			if F <= sphere_f

				holroyd_habit = 's';
				return;
				
			end
			
			% Habit is a graupel
			if F <= graupel_f
			
				holroyd_habit = 'g';
				return;
				
			end
			
			% Habit is a dendrite or plate depending on edge smoothness
			if (F > dendrite_f)
			
				if w_l_m_p < 100 - e_d_r * 200
					holroyd_habit = 'h';
					return;
				else
					holroyd_habit = 'd';
					return;
				end
			
			end
			
			% Habit is an aggregate
			if (F > aggregate_f)
			
				holroyd_habit = 'a';
				return;
				
			end

			% Habit is a plate or a dendrite depending on edge smoothness
			if (((O < 0.3) && (d < 1000)) || (F > plate_f))
				if (d > 700) && (w_l_m_p > 100)
					holroyd_habit = 'd';
					return;
				else
					holroyd_habit = 'h';
					return;
				end
			end
			
			% Habit is irregular (this does include rosettes, which do not seem to be easy to identify with this method)
			holroyd_habit = 'i';
			return;
				
		else
			[x, y, d, a_a, a, r_2, F, O, w_l_p, w_l_m_p, s_c, c_r, e_d_r, m_g, a_r, p] = calc_stat(handles,image_buffer, n_slices);
			
			holroyd_habit = 'C';
			return;
        end
    end
end

%/*************************************************************************/
function [d_length, w_width, d, a_angle, area, r2_correlation, F_fine_detail, O, w_l_r_p, w_l_m_p, sym_coeff, cir_area_ratio, edge_diam_ratio, max_gap, aspect_ratio, p_perimeter_change] = calc_stat(handles, image_buffer, n_slices);

	BITS_PER_SLICE = handles.bits_per_slice;
    MAX_TWOD_DATA_LENGTH = 6000;

	area = 0.0;
	n_count = 0;
	sum_x2= 0.0;
	sum_y2= 0.0;
	sum_x = 0.0;
	sum_y = 0.0;
	sum_xy= 0.0;
	cross_x2= 0.0;
	cross_y2= 0.0;
	cross_xy= 0.0;
	p_perimeter_change = 0;
	min_x = MAX_TWOD_DATA_LENGTH*3;
	min_y = BITS_PER_SLICE*3;
	max_x = 0;
	max_y = 0;
	first_on_edge_x = -1;
	last_on_edge_x = -2;
	first_on_edge_y = -1;
	last_on_edge_y = -2;
	
	d_length = 0;
	w_width = 0;
	a_angle = 0;
	area = 0;
	r2_correlation = 0;
	F_fine_detail = 0;
	O = 0;
	w_l_r_p = 0;
	w_l_m_p = 0;
	sym_coeff = 0;
	cir_area_ratio = 0;
    
	spot_on_off = 0;
	fully_on_count = 0;
	partial_on_count = 0;

	if (n_slices <= 0)
		return;
	end
	
	for i=1:n_slices

		fully_on_temp = 0;
		min_y_temp = BITS_PER_SLICE * 3;
		max_y_temp = 0;
		
		
		for j=1:BITS_PER_SLICE
			if ((image_buffer(i,j)) == '0')     

				tx = i;
				ty =  j;
				if (tx > max_x) 
					max_x = tx;
				end
				
				if (tx < min_x) 
					min_x = tx;
				end
				if (ty > max_y) 
					max_y = ty;
				end
				
				if (ty > max_y_temp)
					max_y_temp = ty;
				end
				
				if (ty < min_y_temp)
					min_y_temp = ty;
				end
				
				if (ty < min_y) 
					min_y = ty;
				end
				
				sum_x2 = sum_x2 + tx * tx;
				sum_y2 = sum_y2 + ty * ty;
				sum_x  = sum_x  + tx;
				sum_y  = sum_y  + ty;
				sum_xy = sum_xy + tx * ty;

                n_count = n_count + 1;
				p(n_count).x = tx;
				p(n_count).y = ty;

				fully_on_temp = fully_on_temp + 1;

				if (spot_on_off == 0)
					spot_on_off = 1;
					p_perimeter_change = p_perimeter_change + 1;
				elseif (j == BITS_PER_SLICE)
					p_perimeter_change = p_perimeter_change + 1;
				end
			else
				if spot_on_off == 1 
					spot_on_off = 0;
					p_perimeter_change = p_perimeter_change + 1;
				end
			end
		end

		if ((image_buffer(i,1)) == '0')
			if first_on_edge_x < 0
				first_on_edge_x = i;
			end
			last_on_edge_x = i;		
		end
			
		if ((image_buffer(i,BITS_PER_SLICE)) == '0')
			if first_on_edge_x < 0
				first_on_edge_x = i;
			end
			
			last_on_edge_x = i;
				
		end

		if (fully_on_temp ~= 0)
			partial_on_count = partial_on_count + 1;
			if (fully_on_temp >= ((max_y_temp - min_y_temp) * 0.8))
				fully_on_count = fully_on_count + 1;
			end
		end
    end
	area = n_count;

%/*** scan the other way for perimeter change ****/

	spot_on_off = 0;
	
	fully_on_count_b = 0;
	partial_on_count_b = 0;
	
	for j=1:BITS_PER_SLICE
	
		fully_on_temp = 0;
		min_x_temp = MAX_TWOD_DATA_LENGTH * 3;
		max_x_temp = 0;
	
		for i=1:n_slices
		
			if ((image_buffer(i,j)) == '0')
			
				tx = i;
		
				if (tx > max_x_temp)
					max_x_temp = ty;
				end
				
				if (tx < min_x_temp)
					min_x_temp = tx;
				end
				
				fully_on_temp = fully_on_temp + 1;
			     
				if (spot_on_off == 0)
					spot_on_off = 1;
					p_perimeter_change = p_perimeter_change + 1;
				elseif (i == n_slices)
					p_perimeter_change = p_perimeter_change + 1;
					
				end
			else
				if (spot_on_off == 1)
					spot_on_off = 0;
					p_perimeter_change = p_perimeter_change + 1;
				end
			end
		end
			
		if ((image_buffer(1,j)) == '0')
			if first_on_edge_y < 0
				first_on_edge_y = j;
			end
			last_on_edge_y = j;
				
		end
		
		if (fully_on_temp ~= 0)
			partial_on_count_b = partial_on_count + 1;
			if (fully_on_temp >= ((max_y_temp - min_y_temp) * 0.8))
				fully_on_count_b = fully_on_count_b + 1;
			end
		end
		
	end

	if (max_x >= min_x) 
		x_length = max_x - min_x +1;
	else
		x_length = 0.0;
	end

	if (max_y >= min_y) 
		y_length = max_y - min_y +1;
	else
		y_length = 0.0;
	end
		
	cross_xy = sum_xy - (sum_x * sum_y / area);
	cross_x2 = sum_x2 - (sum_x * sum_x / area);
	cross_y2 = sum_y2 - (sum_y * sum_y / area);

	slope = cross_xy / cross_x2;
	intercept = (sum_y/(area)) - slope * (sum_x/(area));
	
	r2_correlation = (abs(cross_xy)) / (sqrt( cross_x2 * cross_y2));

	angle_radian = atan(slope);
	a_angle = atan(slope) * (180.0/pi);

	if (a_angle < 0) 
		a_angle = a_angle + 180.0;
		angle_radian = angle_radian + pi;
	end
	
	if (r2_correlation < 0.3)
		a_angle = 0.0;
		angle_radian = 0.0;
	end

	dmin_x = MAX_TWOD_DATA_LENGTH*3;
	dmin_y = BITS_PER_SLICE;
	dmax_x = 0;
	dmax_y = 0;

	if ( (angle_radian > (pi/2.0)) && (angle_radian <= (pi))) 
		angle_radian = (pi - angle_radian);
	elseif ( angle_radian > pi) 
		['HEY: something is wrong here  a_angle = ', num2str(a_angle)]; 
        return
	end
	
	%/*** rotate axes to align with maximum dimension ****/
	
	new_x_vals = [];
	new_y_vals = [];
	
	for i=1:n_count
		new_x = round((p(i).x * cos(angle_radian)) + (p(i).y * sin(angle_radian)));
		new_y = round((p(i).y * cos(angle_radian)) - (p(i).x * sin(angle_radian)));
		
		new_x_vals(end+1) = new_x;
		new_y_vals(end+1) = new_y;
		
		if (new_x > dmax_x) 
			dmax_x = new_x;
		end
		if (new_y > dmax_y) 
			dmax_y = new_y;
		end
		if (new_x < dmin_x) 
			dmin_x = new_x;
		end
		if (new_y < dmin_y) 
			dmin_y = new_y;
		end
	end
	
	new_x_vals_unique = unique(new_x_vals);
	new_y_vals_unique = unique(new_y_vals);
	
	widths_r = [];
	lengths_r = [];
	
	max_gap_x = -32768;
    max_gap_y = -32768;
    prev_x = 32768;
    prev_y = 32768;
    
    %/*** Find particle dimensions for every row and column ****/
	
	for idx=1:length(new_x_vals_unique)
	
		yr_vals_col = [];
	
		curr_x = new_x_vals_unique(idx);
		x_vals_logical = new_x_vals == curr_x;
		
		for idx_idx=1:length(x_vals_logical)
		
			if x_vals_logical(idx_idx) == 1
			
				yr_vals_col(end + 1) = new_y_vals(idx_idx);
				
			end
		
		end
		
		if (length(yr_vals_col)) > 0
        
        	width_r =  max(yr_vals_col) - min(yr_vals_col) + 1;
        
        	widths_r(end+1) = width_r;
        	
        end
        
        if (curr_x - prev_x > max_gap_x)
            max_gap_x = curr_x - prev_x;
        end
        
        prev_x = curr_x;

    end
    
    count_temp = 1;
    
    for idx=1:length(new_y_vals_unique)
	
		xr_vals_row = [];
	
		curr_y = new_y_vals_unique(idx);
		y_vals_logical = new_y_vals == curr_y;
		
		for idx_idx=1:length(y_vals_logical)
		
			if y_vals_logical(idx_idx) == 1
			
				xr_vals_row(end + 1) = new_x_vals(idx_idx);
				
			end
		
		end
        
        if (length(xr_vals_row)) > 0
        
        	length_r =  max(xr_vals_row) - min(xr_vals_row) + 1;
        
        	lengths_r(end+1) = length_r;
        
        end
        
        if (curr_y - prev_y > max_gap_y)
            max_gap_y = curr_y - prev_y;
        end
        
        prev_y = curr_y;

    end
    
    d_length = (dmax_x - dmin_x) + 1;
	w_width = (dmax_y - dmin_y) + 1;
	data_temp = struct('a',{d_length, w_width});
	d = max(data_temp.a);
	
	%/*** Find remaining morphological properties ****/
    
    if (length(lengths_r) > 0) && (length(widths_r) > 0)
    	p_low_x = prctile(lengths_r, 25);
    	p_high_x = prctile(lengths_r, 75);
    	p_higher_x = prctile(lengths_r, 95);
    	var_x = p_high_x - p_low_x;
    	var_higher_x = p_higher_x - p_high_x;
    
    	p_low_y = prctile(widths_r, 25);
    	p_high_y = prctile(widths_r, 75);
    	p_higher_y = prctile(widths_r, 95);
    	var_y = p_high_y - p_low_y;
    	var_higher_y = p_higher_y - p_high_y;
    	
    	y_midpoint = dmax_y - w_width/2;
    	x_midpoint = dmax_x - d_length/2;
    	
    	count_q1 = 0;
    	count_q2 = 0;
    	count_q3 = 0;
    	count_q4 = 0;
    	
    	for idx_idx=1:length(new_x_vals)
    		temp_x_val = new_x_vals(idx_idx);
    		temp_y_val = new_y_vals(idx_idx);
    		
    		if (temp_y_val > y_midpoint) && (temp_x_val > x_midpoint)
    			count_q1 = count_q1 + 1;
    		end
    		if (temp_y_val < y_midpoint) && (temp_x_val > x_midpoint)
    			count_q2 = count_q2 + 1;
    		end
    		if (temp_y_val > y_midpoint) && (temp_x_val < x_midpoint)
    			count_q3 = count_q3 + 1;
    		end
    		if (temp_y_val < y_midpoint) && (temp_x_val < x_midpoint)
    			count_q4 = count_q4 + 1;
    		end
    	end
    	
    	if ((count_q1 == 0) && (count_q4 == 0)) || ((count_q2 == 0) && (count_q3 == 0))
    	
    	    sym_coeff = 0;
    	    
    	else
    	
    		data_temp = struct('a',{count_q1, count_q4});
    	
    		sym_1 = min(data_temp.a)/max(data_temp.a);
    	
    		data_temp = struct('a',{count_q2, count_q3});
    	
    		sym_2 = min(data_temp.a)/max(data_temp.a);
    	
    		sym_coeff = sym_1 * sym_2;
    		
    	end
    
    else
    
    	var_x = 0;
    	var_higher_x = 0;
    	var_y = 0;
    	var_higher_y = 0;
    	sym_coeff = 0;
    
    end
    
    data_temp = struct('a',{max_gap_x, max_gap_y});
    
    max_gap = max(data_temp.a);
    
    max_gap = max_gap - 1;
    
    if max_gap < 0
    	max_gap = 32768;
    end	

    edge_length_y = last_on_edge_y - first_on_edge_y + 1;
    
    edge_length_x = last_on_edge_x - first_on_edge_x + 1;
    
    edge_diam_ratio = edge_length_x/d;
    
    if (edge_length_y > 64)
    	edge_diam_ratio = 1;
    end
    
	data_temp = struct('a',{var_x, var_y});
	
	w_l_r_p = min(data_temp.a);
	
	data_temp = struct('a',{var_higher_x, var_higher_y});
	
	w_l_m_p = min(data_temp.a); 
	
	cir_area_ratio = area/(pi*((d/2)^2));
	
	data_temp = struct('a',{d_length/w_width, w_width/d_length});
    
    if (area < 1.5)
    	aspect_ratio = 1;
    else         
    	aspect_ratio = max(data_temp.a);
    end

	F_fine_detail = p_perimeter_change * d / area;
	
	if (partial_on_count ~=0 ) && (partial_on_count_b ~=0)
		O = ((fully_on_count / partial_on_count) + (fully_on_count_b / partial_on_count_b)) / 2;
	else
		O = 0.0;
	end
end


%/**************************************************************************/
%/*** Not implemented ****/
function result = parabola_fit_center_is_in(image_buffer, n_slices) 


    result = 1;

	return;

end