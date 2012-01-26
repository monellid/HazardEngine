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
		- compute rupture surface
				
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
		#TODO: check lower_seismo_depth > upper_seismo_depth>=0
		# this avoids having a seismogenic layer with thickness equal to 0, and 
		# having a division by zero in case of ruptures with width greater than
		# maximum width allowed by the seismogenic layer.
		#TODO: check minimum(hypocentral depths) >= upper_seismo_depth
		#TODO: check maximum(hypocentral depths) <= lower_seismo_depth
		#TODO: check rupt_aspect_ratio > 0
		self.point = point
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
		mag_index, nodal_plane_index, hypo_depth_index = self.__getMagNodalPlaneHypoDepthIndexes(rupt_index)
		
		occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
		magnitude = occurrence_rates[mag_index][0]
		occurrence_rate = occurrence_rates[mag_index][1] * self.nodal_plane_pmf[nodal_plane_index][3] * self.hypo_depth_pmf[hypo_depth_index][1]
		probability_occurrence = 1 - exp(-occurrence_rate * self.time_span)

		hypocenter = Point(self.point.longitude,self.point.latitude,self.hypo_depth_pmf[hypo_depth_index][0])
		
		strike = self.nodal_plane_pmf[nodal_plane_index][0]
		dip = self.nodal_plane_pmf[nodal_plane_index][1]
		rake = self.nodal_plane_pmf[nodal_plane_index][2]
		
		rup_length, rup_width = self.__getRuptureDimensions(magnitude)

		if rup_length < self.mesh_spacing and rup_width < self.mesh_spacing:
			rup_corners = [hypocenter,hypocenter,hypocenter,hypocenter]
			rupture_surface = numpy.array([hypocenter])
			rupture_surface = rupture_surface.reshape(1,1)
		else:
			# Case 1: rupture width greater than maximum width allowed by
			# seismogenic layer thickness
			maximum_width = (self.lower_seismo_depth - self.upper_seismo_depth) / sin(radians(dip))
			if rup_width > maximum_width:
				rup_length = rup_length * (rup_width / maximum_width)
				rup_width = maximum_width
				w1 = (hypocenter.depth - self.upper_seismo_depth) / sin(radians(dip))
				w2 = (self.lower_seismo_depth - hypocenter.depth) / sin(radians(dip))
			# Case 2: top of rupture depth smaller than upper seismogenic depth
			elif (rup_width / 2) > ((hypocenter.depth - self.upper_seismo_depth) / sin(radians(dip))):
				w1 = (hypocenter.depth - self.upper_seismo_depth) / sin(radians(dip))
				w2 = rup_width - w1
			# Case 3: bottom of rupture depth greater then lower seismogenic depth
			elif (rup_width / 2) > ((self.lower_seismo_depth - hypocenter.depth) / sin(radians(dip))):
				w2 = (self.lower_seismo_depth - hypocenter.depth) / sin(radians(dip))
				w1 = rup_width - w2
			# Case 4: rupture inside seismogenic layer
			else: 
				w1 = (rup_width / 2)
				w2 = (rup_width / 2)

			upper_hor_distance = sqrt( (w1 * cos(radians(dip)))**2 + (rup_length / 2)**2 )
			upper_vert_distance = - w1 * sin(radians(dip))
			upper_right_azimuth = strike - degrees(atan( (2 * w1 * cos(radians(dip))) / rup_length ))
			upper_left_azimuth = strike + 180.0 + degrees(atan( (2 * w1 * cos(radians(dip))) / rup_length ))

			upper_right = hypocenter.getPoint(upper_hor_distance,upper_vert_distance,upper_right_azimuth)
			upper_left = hypocenter.getPoint(upper_hor_distance,upper_vert_distance,upper_left_azimuth)

			lower_hor_distance = sqrt( (w2 * cos(radians(dip)))**2 + (rup_length / 2)**2 )
			lower_vert_distance = + w2 * sin(radians(dip))
			lower_right_azimuth = strike + degrees(atan( (2 * w2 * cos(radians(dip))) / rup_length ))
			lower_left_azimuth = strike + 180.0 - degrees(atan( (2 * w2 * cos(radians(dip))) / rup_length ))

			lower_right = hypocenter.getPoint(lower_hor_distance,lower_vert_distance,lower_right_azimuth)
			lower_left = hypocenter.getPoint(lower_hor_distance,lower_vert_distance,lower_left_azimuth)

			rup_corners = [upper_left,upper_right,lower_left,lower_right]
			rupture_surface = self.__getRuptureSurface(rup_corners)

		return {'magnitude':magnitude,'strike':strike,'dip':dip,'rake':rake,
				'tectonic':self.tectonic_region_type,'hypocenter':hypocenter,
				'upper_left':rup_corners[0],'upper_right':rup_corners[1],
				'lower_left':rup_corners[2],'lower_right':rup_corners[3],
				'surface':rupture_surface,
				'rate':occurrence_rate,'probability':probability_occurrence}
			
	def __getMagNodalPlaneHypoDepthIndexes(self,rupt_index):
		num_mags = len(self.freq_mag_dist.getAnnualOccurrenceRates())
		num_nodal_planes = len(self.nodal_plane_pmf)
		num_hypo_depths = len(self.hypo_depth_pmf)
		
		mag_index = rupt_index % num_mags
		nodal_plane_index = (rupt_index / num_mags) % num_nodal_planes
		hypo_depth_index = (rupt_index / (num_mags * num_nodal_planes)) % num_hypo_depths
		
		return mag_index, nodal_plane_index, hypo_depth_index
			
	def __getRuptureDimensions(self,magnitude):
		area = self.mag_scaling_rel.getMedianArea(magnitude)
		rup_length = sqrt(area * self.rupt_aspect_ratio)
		rup_width = area / rup_length
		return rup_length, rup_width
			
	def __getRuptureSurface(self,rup_corners):
		"""
		Compute surface mesh from rupture coordinates.
		"""
		top_edge = rup_corners[0].getEquallySpacedPoints(rup_corners[1],self.mesh_spacing)
		bottom_edge = rup_corners[2].getEquallySpacedPoints(rup_corners[3],self.mesh_spacing)
		mesh_points = []
		for i in range(len(top_edge)):
			top = top_edge[i]
			bottom = bottom_edge[i]
			points = top.getEquallySpacedPoints(bottom,self.mesh_spacing)
			mesh_points.extend(points)
		# organize mesh points into a 2D array
		# number of rows corresponds to number of points along dip
		# number of columns corresponds to number of points along strike
		surface = numpy.array(mesh_points)
		surface = surface.reshape(len(top_edge),len(mesh_points)/len(top_edge))
		surface = numpy.transpose(surface)
		
		return surface