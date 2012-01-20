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
		self.tectonic_region_type = tectonic_region_type
		self.time_span = time_span
		self.mesh_spacing = mesh_spacing
		
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
			mag_index, nodal_plane_index, hypo_depth_index = __getMagNodalPlaneHypoDepthIndexes(self,rup_index)
			
			occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
			magnitude = occurrence_rates[mag_index][0]
			occurrence_rate = occurrence_rates[mag_index][1] * self.nodal_plane_pmf[nodal_plane_index][3] * self.hypo_depth_pmf[hypo_depth_index][1]
			probability_occurrence = 1 - exp(-occurrence_rate * time_span)

			hypocenter = Point(self.point.longitude,self.point.latitude,self.hypo_depth_pmf[hypo_depth_index][0])
			
			strike = self.nodal_plane_pmf[nodal_plane_index][0]
			dip = self.nodal_plane_pmf[nodal_plane_index][1]
			rake = self.nodal_plane_pmf[nodal_plane_index][2]
			
			rup_length, rup_width = __getRuptureDimensions(self,magnitude)
			
			if rup_length < self.mesh_spacing and rup_width < self.mesh_spacing:
				rupture_surface = numpy.array([hypocenter])
			else:
				rup_corners = self.__getRuptureCorners(hypocenter,rup_length,rup_width,strike,dip)
				rupture_surface = __getRuptureSurface(self,rup_corners)
				
			return {'magnitude':magnitude,'strike':strike,'dip':dip,'rake':rake,
					'hypocenter':hypocenter,'surface':rupture_surface,
					'rate':occurrence_rate,'probability':probability_occurrence}
				
		def __getMagNodalPlaneHypoDepthIndexes(self,rup_index):
			num_mags = len(self.freq_mag_dist.getAnnualOccurrenceRates())
			num_nodal_planes = len(self.nodal_plane_pmf)
			num_hypo_depth = len(self.hypo_depth_pmf)
			
			mag_index = rupt_index % num_mag_values
			nodal_plane_index = (rupt_index / num_mag_values) % num_nodal_planes
			hypo_depth_index = (rupt_index / (num_mag_values * num_nodal_planes)) % num_hypo_depths
			
			return mag_index, nodal_plane_index, hypo_depth_index
				
		def __getRuptureDimensions(self,magnitude):
			area = self.mag_scaling_rel.getMedianArea(magnitude)
			rup_length = sqrt(area * self.rup_aspect_ratio)
			rup_width = area / rup_length
			return rup_length, rup_width
				
		def __getRuptureSurface(self,rup_corners):
			"""
			Compute surface mesh from rupture coordinates.
			"""
			top_edge = rup_corners[0].getEquallySpacedPoints(rup_corners[1],self.mesh_spacing)
			bottom_edge = rup_corners[2].getEquallySpacedPoints(rup_corners[3],self.mesh_spacing)
			mesh_points = []
			for top,bottom in top_edge,bottom_edge:
				points = top.getEquallySpacedPoints(bottom,self.mesh_spacing)
				mesh_points.extend(points)
			# organize mesh points into a 2D array
			# number of rows corresponds to number of points along dip
			# number of columns corresponds to number of points along strike
			surface = numpy.array(mesh_points)
			surface = surface.reshape(len(mesh_points)/len(top_edge),len(top_edge))
			
			return surface
			
				
		def __getRuptureCorners(self,hypocenter,rup_length,rup_width,strike,dip):
			"""
			From rupture's hypocenter, length, width, strike and dip computes rupture's corner coordinates
			(assuming rupture hypocenter being the rupture surface centroid)
			"""
			strike = radians(strike)
			dip = radians(dip)
			
			aspect_ratio = rup_length / rup_width
			horizontal_distance = sqrt( (rup_width * cos(dip) / 2)**2 + (rup_length / 2)**2 )
		 	vert_distance = rup_width * sin(dip) / 2
			# left and right are with respect to rupture hypocenter
			azimuth_right = degrees(strike - atan( cos(dip) / aspect_ratio ))
			azimuth_left = degrees(azimuth_right + 180 + 2 * atan( cos(dip) / aspect_ratio ))
			
			upper_left = hypocenter.getPoint(horizontal_distance,-vert_distance,azimuth_left)
			upper_right = hypocenter.getPoint(horizontal_distance,-vert_distance,azimuth_right)
			lower_left = hypocenter.getPoint(horizontal_distance,+vert_distance,azimuth_left)
			lower_right = hypocenter.getPoint(horizontal_distance,+vert_distance,azimuth_right)
			
			return self.__adjustRuptureCorners([upper_left,upper_right,lower_left,lower_right],strike,dip)
			
		def __adjustRuptureCorners(self,rup_corners,strike,dip):
			"""
			Adjusts rupture's corners so that rupture fits inside the seismogenic layer.
			"""
			# If the rupture's top edge depth is smaller than the upper seismogenic depth or
			# the rupture's bottom edge depth is larger than the lower seismogenic depth,
			# shifts the rupture (along the dip direction) downwards or upwards respectively.
			# if the resulting bottom or top edges goes off of the seismogenic layer, the
			# rupture is reshaped: rupture width is clipped to the width allowed by the
			# seismogenic layer thickness, and extended symmetrically along the strike
			# direction. Reshaping is done so as rupture area is conserved.
			if rup_corners[0].depth < self.upper_seismo_depth:
				vertical_increment = self.upper_seismo_depth - rup_corners[0].depth
				rup_corners = __shiftRuptureCornersAlongDipDirection(rup_corners,vertical_increment,strike,dip)
				if rup_corners[2].depth > self.lower_seismo_depth:
					delta_width = (rup_corners[2].depth - lower_seismo_depth) / sin(dip)
					delta_length = (rup_length * delta_width) / (2 * rup_width)
					vertical_dist = self.lower_seismo_depth - self.upper_seismo_depth
					horizontal_dist = vertical_dist / tan(dip)
					rup_corners[0] = rup_corners[0].getPoint(delta_length,0.0,strike + 180.0)
					rup_corners[1] = rup_corners[1].getPoint(delta_length,0.0,strike)
					rup_corners[2] = rup_corners[0].getPoint(horizontal_dist,vertical_dist,strike + 90.0)
					rup_corners[3] = rup_corners[1].getPoint(horizontal_dist,vertical_dist,strike + 90.0)
			if rup_corners[2].depth > self.lower_seismo_depth:
				vertical_increment =  self.lower_seismo_depth - rup_corners[2].depth
				rup_corners =  __shiftRuptureCornersAlongDipDirection(rup_corners,vertical_increment,strike,dip)
				if rup_corners[0].depth < upper_seismo_depth:
					delta_width = (upper_seismo_depth - rup_corners[0].depth) / sin(dip)
					delta_length = (rup_length * delta_width) / (2 * rup_width)
					vertical_dist = self.lower_seismo_depth - self.upper_seismo_depth
					horizontal_dist = vertical_dist / tan(dip)
					rup_corners[2] = rup_corners[2].getPoint(delta_length,0.0,strike + 180.0)
					rup_corners[3] = rup_corners[3].getPoint(delta_length,0.0,strike)
					rup_corners[0] = rup_corners[2].getPoint(horizontal_dist,-vertical_dist,strike + 270.0)
					rup_corners[1] = rup_corners[3].getPoint(horizontal_dist,-vertical_dist,strike + 270.0)
					
		def __shiftRuptureCornersAlongDipDirection(rup_corners,vertical_increment,strike,dip):
			"""
			The method assumes the 4 rupture corners to be aligned on a plane with an
			inclination with respect to the Earth surface equal to dip, and the top and
			bottom corners to be aligned along the strike direction. That is the 4
			corners must describe a planar rectangle.
			If vertical_increment is positive, corners are shifted downwards.
			If vertical_increment is negative, corners are shifted upwards.
			"""
			if vertical_increment > 0:
				azimuth = strike + 90.0
			else:
				azimuth = strike + 90.0 + 180.0
			horizontal_distance = abs(vertical_increment) / tan(dip)
			return [point.getPoint(horizontal_distance,vertical_increment,azimuth) for point in rup_corners]