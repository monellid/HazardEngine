from geo import *

class PoissonianPointSource:
	"""
	Class defining Poissonian point source. Ruptures are centered on
	a single geographical location and have Poissonian probability of
	occurrence.
	"""
	
	def __init__(self,point,freq_mag_dist,nodal_plane_pmf,hypo_depth_pmf,upper_seismo_depth,lower_seismo_depth,mag_scaling_rel,rupt_aspect_ratio,tectonic_region_type,time_span,mesh_spacing):
		"""
		point: geographical location at the Earth surface corresponding to surface projection of ruptures' centroids.
		freq_mag_dist: frequency magnitude distribution
		nodal_plane_dist: nodal plane (strike,dip,rake) probability mass function [(strike1,dip1,rake1,weight1),...,(striken,dipn,raken,weightn)]
		hypo_depth_pmf: hypocentral depth (in km) probability mass function [(hypo1,weight1),...,(hypo2,weight2)]
		upper_seismo_depth: upper seismogenic depth (in km)
		lower_seismo_depth: lower seismogenic depth (in km)
		mag_scaling_rel: magnitude-area scaling relationship
		rupt_aspect_ratio: rupture aspect ratio
		tectonic_region_type: tectonic region type
		time_span: time span (years) for calculating probability of occurrence
		mesh_spacing: rupture surface mesh spacig (in km)
		
		Ruptures are generated as follows:
		Loop over magnitude values in the freq_mag_dist. For each mag value,
		loop over nodal planes. For each nodal plane, loop over hypocentral
		depths. For each depth define a rupture in the following way:
		- rupture_hypocenter = Point(point.longitude,point.latitude,current_depth)
		- magnitude = current_magnitude_value
		- strike = current_strike
		- dip = current_dip
		- rake = current_rake
		- occurrence_rate = current_occurrence_rate * nodal_plane_probability * hypocentral_depth_probability
		- probability_occurrence = 1 - exp(-occurrence_rate * time_span)
		- area = mag_scaling_rel.getMedianArea(magnitude)
		- rup_length = sqrt(area * rup_aspect_ratio)
		- rup_width = area / rup_length
		- compute rupture surface:
			if rup_length < mesh_spacing and rup_width < mesh_spacing:
				rupture_surface = rupture_hypocenter
			else
				
				From the rupture hypocenter, computes rupture's corners coordinates as follows:
				- horizontal_distance = sqrt( (rup_width * cos(dip) / 2)**2 + (rup_length / 2)**2 )
				- vert_increment = rup_width * sin(dip) / 2
				- azimuth_right = strike - atan( cos(dip) / aspect_ratio )
				- azimuth_left = azimuth_right + 180 + 2 * atan( cos(dip) / aspect_ratio )
				# left and right are with respect to rupture hypocenter
				- rupture_upper_left_corner = rupture_hypocenter.getPoint(horizontal_distance,-vertical_increment,azimuth_left)
				- rupture_upper_right_corner = rupture_hypocenter.getPoint(horizontal_distance,-vertical_increment,azimuth_right)
				- rupture_lower_left_corner = rupture_hypocenter.getPoint(horizontal_distance,+vertical_increment,azimuth_left)
				- rupture_lower_right_corner = rupture_hypocenter.getPoint(horizontal_distance,+vertical_increment,azimuth_right)
				
				# Check if rupture fits inside seismogenic layer:
				if rupture_upper_left_corner.depth < upper_seismo_depth:
					# shift downwards until top edge reaches the upper seismogenic depth,
					# and if the resulting bottom edge goes deeper than lower seismogenic depth, reshape the rupture
					azimuth = rupture_upper_left_corner.getAzimuth(rupture_lower_left_corner)
					vertical_increment = upper_seismo_depth - rupture_upper_left_corner.depth
					horizontal_distance = vertical_increment / tan(dip)
					rupture_upper_left_corner = rupture_upper_left_corner.getPoint(horizontal_distance,vertical_increment,azimuth)
					rupture_upper_right_corner = rupture_upper_right_corner.getPoint(horizontal_distance,vertical_increment,azimuth)
					rupture_lower_left_corner = rupture_lower_left_corner.getPoint(horizontal_distance,vertical_increment,azimuth)
					rupture_lower_right_corner = rupture_lower_right_corner.getPoint(horizontal_distance,vertical_increment,azimuth)
					if rupture_lower_left_corner.depth > lower_seismo_depth:
						# clip fault plane up to the lower seismogenic depth, and extends it along length
						# compute width offset
						delta_width = (rupture_lower_left_corner.depth - lower_seismo_depth) / sin(dip)
						# compute length offset
						delta_length = (rup_length * delta_width) / (2 * rup_width)
						# recompute new upper corners but keep the azimuths first
						azimuth_along_length_clockwise = rupture_upper_left_corner.getAzimuth(rupture_upper_right_corner)
						azimuth_along_width_downwards = rupture_upper_left_corner.getAzimuth(rupture_lower_left_corner)
						rupture_upper_left_corner = rupture_upper_left_corner.getPoint(delta_length,0.0,azimuth_along_length_clockwise + 180.0)
						rupture_upper_right_corner = rupture_upper_right_corner.getPoint(delta_length,0.0,azimuth_along_length_clockwise)
						vertical_dist = lower_seismo_depth - upper_seismo_depth
						horizontal_dist = vertical_dist / tan(dip)
						rupture_lower_left_corner = rupture_upper_left_corner.getPoint(horizontal_dist,vertical_dist,azimuth_along_width_downwards)
						rupture_lower_right_corner = rupture_upper_right_corner.getPoint(horizontal_dist,vertical_dist,azimuth_along_width_downwards)
						
					
				if rupture_lower_left_corner.depth > lower_seismo_depth:
					# shift upwards until bottom edge reaches the lower seismogenic depth,
					# and if the resulting top edge goes shallower then upper seismogenic depth, reshape the rupture
					azimuth = rupture_lower_left_corner.getAzimuth(rupture_upper_left_corner)
					vertical_distance =  rupture_lower_left_corner.depth - lower_seismo_depth
					horizontal_distance = vertical_distance / tan(dip)
					rupture_upper_left_corner = rupture_upper_left_corner.getPoint(horizontal_distance,-vertical_distance,azimuth)
					rupture_upper_right_corner = rupture_upper_right_corner.getPoint(horizontal_distance,-vertical_distance,azimuth)
					rupture_lower_left_corner = rupture_lower_left_corner.getPoint(horizontal_distance,-vertical_distance,azimuth)
					rupture_lower_right_corner = rupture_lower_right_corner.getPoint(horizontal_distance,-vertical_distance,azimuth)
					if rupture_upper_left_corner.depth < upper_seismo_depth:
						# clip fault plane up to the lower seismogenic depth, and extends it along length
						# compute width offset
						delta_width = (upper_seismo_depth - rupture_upper_left_corner.depth) / sin(dip)
						# compute length offset
						delta_length = (rup_length * delta_width) / (2 * rup_width)
						# recompute new upper corners but keep the azimuths first
						azimuth_along_length_clockwise = rupture_upper_left_corner.getAzimuth(rupture_upper_right_corner)
						azimuth_along_width_upwards = rupture_lower_left_corner.getAzimuth(rupture_upper_left_corner)
						rupture_lower_left_corner = rupture_lower_left_corner.getPoint(delta_length,0.0,azimuth_along_length_clockwise + 180.0)
						rupture_lower_right_corner = rupture_lower_right_corner.getPoint(delta_length,0.0,azimuth_along_length_clockwise)
						vertical_dist = lower_seismo_depth - upper_seismo_depth
						horizontal_dist = vertical_dist / tan(dip)
						rupture_upper_left_corner = rupture_lower_left_corner.getPoint(horizontal_dist,-vertical_dist,azimuth_along_width_upwards)
						rupture_upper_right_corner = rupture_lower_right_corner.getPoint(horizontal_dist,-vertical_dist,azimuth_along_width_upwards)
				
				To construct the rupture surface, compute rupture top edge coordinates (by resampling 
				[using method getEquallySpacedPoints] the line connecting rupture_upper_left_corner 
				with rupture_upper_right_corner with a distance equal to mesh_spacing). For each point on the rupture
				top edge compute correspondig point on the rupture bottom edge,that is:
				mesh_points = []
				for point in rupture_top_edge:
					vertical_increment = rup_width * sin(dip)
					horizontal_distance = vertical_increment / tan(dip) NOTE: for dip=90 put horizontal_distance = 0 without the formula
					azimuth = strike + 90.0
					bottom_edge_point = top_edge_point.getPoint(horizontal_distance,vertical_increment,azimuth)
					points = point.getEquallySpacedPoints(bottom_edge_point,mesh_spacing)
					mesh_points.extend(points)
				# organize mesh points into a 2D array
				# number of rows corresponds to number of points along dip
				# number of columns corresponds to number of points along strike
				surface = numpy.array(mesh_points)
				surface = surface.reshape(len(top_edge),len(mesh_points)/len(top_edge))
				surface = numpy.transpose(surface)
				
		The total number of ruptures can be computed as follows:
		- total_num_rupture = num_mag_values * num_nodal_planes * num_hypo_depths
		
		To get rupture of index i:
		mag_index = i % num_mag_values
		nodal_plane_index = (i / num_mag_values) % num_nodal_planes
		hypo_depth_index = (i / (num_mag_values * num_nodal_planes)) % num_hypo_depths
		"""
		#TODO: check point is at the Earth surface (depth==0)
		#TODO: check PMFs. probability values must be positive and sum up to 1
		#TODO: check strike, dip, and rake values (0<=strike<=360, 0<dip<=90,-180<=rake<=180)
		#TODO: check hypocentral depths. All depths must be positive and greater than 0.
		#TODO: check lower_seismo_depth>=upper_seismo_depth>=0
		#TODO: check minimum(hypocentral depths) >= upper_seismo_depth
		#TODO: check maximum(hypocentral depths) <= lower_seismo_depth
		#TODO: check rupt_aspect_ratio > 0
		self.freq_mag_dist = freq_mag_dist
		self.nodal_plane_pmf = nodal_plane_pmf
		self.hypo_depth_pmf = hypo_depth_pmf
		self.upper_seismo_depth = upper_seismo_depth
		self.lower_seismo_depth = lower_seismo_depth
		self.mag_scaling_rel = mag_scaling_rel
		self.rupt_aspect_ratio = rupt_aspect_ratio
		
		def getNumRuptures(self):
			"""
			Return number of ruptures
			"""
			num_mags = len(self.freq_mag_dist.getAnnualOccurrenceRates())
			num_nodal_planes = len(self.nodal_plane_pmf)
			num_hypo_depth = len(self.hypo_depth_pmf)
			
			return num_mags * num_nodal_planes * num_hypo_depth
	
		def getRupture(self,rupt_index):
			"""
			Returns rupture for
			index rupt_index.
			"""
			occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
			num_mags = len(occurrence_rates)
			num_nodal_planes = len(self.nodal_plane_pmf)
			num_hypo_depth = len(self.hypo_depth_pmf)
			
			mag_index = rupt_index % num_mag_values
			nodal_plane_index = (rupt_index / num_mag_values) % num_nodal_planes
			hypo_depth_index = (rupt_index / (num_mag_values * num_nodal_planes)) % num_hypo_depths
			
			rupture_hypocenter = Point(self.point.longitude,self.point.latitude,self.hypo_depth_pmf[hypo_depth_index][0])
			magnitude = occurrence_rates[mag_index][0]
			strike = self.nodal_plane_pmf[nodal_plane_index][0]
			dip = self.nodal_plane_pmf[nodal_plane_index][1]
			rake = self.nodal_plane_pmf[nodal_plane_index][2]
			occurrence_rate = occurrence_rates[mag_index][1] * self.nodal_plane_pmf[nodal_plane_index][3] * self.hypo_depth_pmf[hypo_depth_index][1]
			probability_occurrence = 1 - exp(-occurrence_rate * time_span)
			area = mag_scaling_rel.getMedianArea(magnitude)
			rup_length = sqrt(area * rup_aspect_ratio)
			rup_width = area / rup_length