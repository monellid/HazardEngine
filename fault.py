from math import *
from geo import *
from rup import *
import numpy
		
class SimpleFaultSurface:
	"""
	Class defining Simple Fault Surface.
	"""
	
	def __init__(self,fault_trace,upper_seismo_depth,lower_seismo_depth,dip,mesh_spacing = 1):
		"""
		Represents fault surface as regular (uniformly spaced) 3D mesh of points derived from:
		fault_trace: list of points representing intersection between fault surface and Earth surface
		upper_seismo_depth: upper seismogenic depth (i.e. depth to fault's top edge, in km)
		lower_seismo_depth: lower seismogenic depth (i.e. depth to fault's bottom edge, in km)
		dip: dip angle (in degrees)
		mesh_spacing: spacing (in km) between mesh points
		"""
		# TODO: check that all points in fault_trace have depth = 0
		# TODO: check that there are at least two points on the fault trace and they are not coincidents
		# TODO: check upper_seismo_depth >= 0
		# TODO: check (lower_seismo_depth - upper_seismo_depth) / sin(dip) >= mesh_spacing
		# TODO: check dip >0 && dip <=90
		# TODO: check mesh_spacing > 0
		self.fault_trace = fault_trace
		self.upper_seismo_depth = upper_seismo_depth
		self.lower_seismo_depth = lower_seismo_depth
		self.dip = dip
		self.mesh_spacing = mesh_spacing
		self.surface = self.getSurfaceMesh()
		
	def getSurfaceMesh(self):
		"""
		Computes 3D mesh representing fault surface. The mesh is constructed by 
		translating the fault trace along the dip angle, along a direction perpendicular
		to the fault trace strike (computed considering the first and last points in the fault trace).
		Distance between mesh points is given by the mesh_spacing parameter (in km).
		"""		
		top_edge = self.__getFaultTopEdge()
		
		# compute mesh points. Loops over points in the top edge, for each point
		# on the top edge compute corresponding point on the bottom edge, then
		# computes equally spaced points between top and bottom points
		vertical_distance = (self.lower_seismo_depth - self.upper_seismo_depth)
		horizontal_distance = vertical_distance / tan(radians(self.dip))
		strike = self.fault_trace[0].getAzimuth(self.fault_trace[-1])
		azimuth = strike + 90.0
		mesh_points = []
		for point in top_edge:
			bottom_point = point.getPoint(horizontal_distance,vertical_distance,azimuth)
			points = point.getEquallySpacedPoints(bottom_point,self.mesh_spacing)
			mesh_points.extend(points)
		
		# organize mesh points into a 2D array
		# number of rows corresponds to number of points along dip
		# number of columns corresponds to number of points along strike
		surface = numpy.array(mesh_points)
		surface = surface.reshape(len(top_edge),len(mesh_points)/len(top_edge))
		surface = numpy.transpose(surface)
		return surface
		
	def __getFaultTopEdge(self):
		"""
		Computes fault top edge coordinates by translating the fault trace along dip
		from the Earth surface to the upper seismogenic depth in a direction 
		perpendicular to the fault strike (computed considering the first and last points).
		Point's coordinates are then resampled to be equally spaced (with spacing equal
		to mesh_spacing).
		"""
		top_edge = []
		
		if self.dip < 90.0:
			horizontal_distance = self.upper_seismo_depth / tan(radians(self.dip))
		else:
			horizontal_distance = 0.0
		vertical_distance = self.upper_seismo_depth
		strike = self.fault_trace[0].getAzimuth(self.fault_trace[-1])
		azimuth = strike + 90.0
		for point in self.fault_trace:
			top_edge.append(point.getPoint(horizontal_distance,vertical_distance,azimuth))

		top_edge = Line(top_edge).getResampledLine(self.mesh_spacing).point_list
		return top_edge
		
class PoissonianFaultSource:
	"""
	Class defining Poissonian Fault source, that is a fault source with
	ruptures having Poissonian probability of occurrence.
	"""
	
	def __init__(self,fault_surf,freq_mag_dist,mag_scaling_rel,rake,rup_aspect_ratio,tectonic_region_type,time_span):
		"""
		fault_surf: fault surface
		freq_mag_dist: frequency magnitude distribution
		mag_scaling_rel: magnitude scaling relationship
		rake: rake angle (-180 <= rake <= 180)
		rup_aspect_ratio: rupture aspect ratio (> 0)
		tectonic_region_type: tectonic region type
		time_span: time span (>= 0)
		"""
		#TODO: check rake >= -180 & rake <=180
		#TODO: check rup_aspect_ratio > 0
		#TODO: check time_span >= 0
		self.fault_surf = fault_surf
		self.freq_mag_dist = freq_mag_dist
		self.mag_scaling_rel = mag_scaling_rel
		self.rake = rake
		self.rup_aspect_ratio = rup_aspect_ratio
		self.tectonic_region_type = tectonic_region_type
		self.time_span = time_span
		self.rupture_data = self.__getRuptureData()
		
	def __getRuptureData(self):
		"""
		Return list containing rupture data. The length of the list
		corresponds to the number of ruptures. Each entry in the list
		is a dictionary with the following keys:
		- 'mag': rupture magnitude
		- 'rate': rupture annual occurrence rate
		- 'first': tuple (i,j) containing indexes of rupture's first mesh point
		- 'last_length': tuple (i,j) containing indexes of rupture's last mesh point along length
		- 'last_width': tuple (i,j) containing indexes of rupture's last mesh point along width
		"""
		rupture_data = []
		
		# Loop over magnitude values in the frequency magnitude distribution.
		# For each magnitude, computes rupture's dimensions (length and width).
		# Based on rupture's dimensions and rupture's offset (assumed equal
		# to surface mesh spacing), computes rupture's first and last point
		# indexes.
		num_points_along_length = self.fault_surf.surface.shape[1]
		num_points_along_width = self.fault_surf.surface.shape[0]
		fault_length = (num_points_along_length - 1) * self.fault_surf.mesh_spacing
		fault_width = (num_points_along_width - 1) * self.fault_surf.mesh_spacing
		occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
		for mag,rate in occurrence_rates:
			
			# compute rupture dimensions
			area = self.mag_scaling_rel.getMedianArea(mag)
			rup_length = sqrt(area * self.rup_aspect_ratio)
			rup_width = area / rup_length
			
			# reshape rupture (conserving area) if its length or width
			# exceeds fault's lenght or width
			if rup_length < fault_length and rup_width > fault_width:
				rup_length = rup_length * (rup_width / fault_width)
				rup_width = fault_width
			elif rup_length > fault_length and rup_width < fault_width:
				rup_width = rup_width * (rup_length / fault_length)
				rup_length = fault_length
				
			# clip rupture's length and width to
			# fault's length and width if both rupture
			# dimensions are greater or equal than fault dimensions
			if rup_length >= fault_length and rup_width >= fault_width:
				rup_length = fault_length
				rup_width = fault_width
			
			# round rupture dimensions with respect to mesh_spacing and compute number
			# of points in the rupture along length and strike
			rup_length = round(rup_length / self.fault_surf.mesh_spacing) * self.fault_surf.mesh_spacing
			rup_width = round(rup_width / self.fault_surf.mesh_spacing) * self.fault_surf.mesh_spacing
			num_rup_points_along_length = int(rup_length / self.fault_surf.mesh_spacing) + 1
			num_rup_points_along_width = int(rup_width / self.fault_surf.mesh_spacing) + 1
							
			# based on rupture dimensions, compute number of ruptures along length and width
			# for each rupture compute first and last points
			if rup_length < self.fault_surf.mesh_spacing and rup_width < self.fault_surf.mesh_spacing:
				num_rup_along_length = num_points_along_length
				num_rup_along_width = num_points_along_width					
			elif rup_length == fault_length and rup_width == fault_width:
				num_rup_along_length = 1
				num_rup_along_width = 1
			else:
				num_rup_along_length = int(int(fault_length - rup_length) / self.fault_surf.mesh_spacing + 1)
				num_rup_along_width = int(int(fault_width - rup_width) / self.fault_surf.mesh_spacing + 1)
			num_rup = num_rup_along_length * num_rup_along_width
			for i in range(num_rup_along_width):
				for j in range(num_rup_along_length):
					first_point = (i,j)
					last_point_along_length = (i,j + num_rup_points_along_length - 1)
					last_point_along_width = (i + num_rup_points_along_width - 1, j)
					data = {'mag':mag,
							'rate':rate / num_rup,
							'first':first_point,
							'last_length':last_point_along_length,
							'last_width':last_point_along_width}
					rupture_data.append(data)
		return rupture_data
					
	def getNumRuptures(self):
		"""
		Return number of ruptures.
		"""
		return len(self.rupture_data)
		
	def getRupture(self,rupt_index):
		"""
		Return rupture corresponding to rupt_index.
		A rupture is currently defined in terms of:
		- magnitude
		- strike (defined as the azimuth between the first and last locations of the rupture top edge for extended ruptures
					 or as azimuth between first and last locations of fault trace for point ruptures)
		- dip (as derived from the fault surface)
		- rake
		- tectonic region type
		- hypocenter (defined as the centroid of the rupture surface)
		- rupture surface mesh
		- rate of occurrence
		- probability of occurrence
		"""
		mag = self.rupture_data[rupt_index]['mag']
		rate = self.rupture_data[rupt_index]['rate']
		
		first = self.rupture_data[rupt_index]['first']
		last_length = self.rupture_data[rupt_index]['last_length']
		last_width = self.rupture_data[rupt_index]['last_width']
		rup_surf_mesh = self.fault_surf.surface[first[0]:last_width[0]+1,first[1]:last_length[1]+1]
		
		# if point source measure strike as azimuth from first and last points in the fault trace
		if len(rup_surf_mesh)==1:
			strike = self.fault_surf.fault_trace[0].getAzimuth(self.fault_surf.fault_trace[-1])
		else:
			# compute strike as azimuth from first and last points in rupture top edge
			strike = rup_surf_mesh[0,0].getAzimuth(rup_surf_mesh[0,-1])

		hypocenter = self.__getSurfaceCentroid(rup_surf_mesh)
		
		# Poissonian probability of one or more occurrences
		probability_occurrence = 1 - exp(-rate * self.time_span)
		
		return {'magnitude':mag,'strike':strike,'dip':self.fault_surf.dip,'rake':self.rake,
				'tectonic':self.tectonic_region_type,'hypocenter':hypocenter,
				'surface':rup_surf_mesh,
				'rate':rate,'probability':probability_occurrence}
		
	def __getSurfaceCentroid(self,surf):
		"""
		Computes surface centroid, that is its geometrical center.
		"""
		# if number of grid points along length and width is odd, returns
		# location corresponding to central index
		# if number of grid points along length is even and along width is
		# odd (or viceversa), returns centroid of the central segment.
		# if number of grid points along length and width is even, returns
		# centroid of the central patch
		num_rows = surf.shape[0]
		num_cols = surf.shape[1]
		if isOdd(num_rows) is True and isOdd(num_cols) is True:
			return surf[int((num_rows - 1) / 2),int((num_cols - 1) / 2)]
		elif isOdd(num_rows) is True and isOdd(num_cols) is not True:
			p1 = surf[int((num_rows - 1) / 2),int(num_cols / 2 - 1)]
			p2 = surf[int((num_rows - 1) / 2),int(num_cols / 2)]
			mean_lon = (p1.longitude + p2.longitude) / 2
			mean_lat = (p1.latitude + p2.latitude) / 2
			mean_depth = (p1.depth + p2.depth) / 2
			return Point(mean_lon,mean_lat,mean_depth)
		elif isOdd(num_rows) is not True and isOdd(num_cols) is True:
			p1 = surf[int(num_rows / 2 - 1),int((num_cols - 1) / 2)]
			p2 = surf[int(num_rows / 2),int((num_cols - 1)) / 2]
			mean_lon = (p1.longitude + p2.longitude) / 2
			mean_lat = (p1.latitude + p2.latitude) / 2
			mean_depth = (p1.depth + p2.depth) / 2
			return Point(mean_lon,mean_lat,mean_depth)
		else:
			p1 = surf[int(num_rows / 2 - 1),int(num_cols / 2 - 1)]
			p2 = surf[int(num_rows / 2 - 1),int(num_cols / 2)]
			p3 = surf[int(num_rows / 2),int(num_cols / 2)]
			p4 = surf[int(num_rows / 2),int(num_cols / 2 - 1)]
			mean_lon = (p1.longitude + p2.longitude + p3.longitude + p4.longitude) / 4
			mean_lat = (p1.latitude + p2.latitude + p3.latitude + p4.latitude) / 4
			mean_depth = (p1.depth + p2.depth + p3.depth + p4.depth) / 4
			return Point(mean_lon,mean_lat,mean_depth)
		
	def getMinimumDistance(self,point):
		"""
		Compute minimum (horizontal) distance between source and given point.
		The minimum distance is defined as the minimum great circle
		distance between the given point and the surface projection of the
		points constituting the fault surface boundary.
		"""
		# TODO: if fault is vertical (dip = 90) compute distance only to
		# surface projections of points on the top edge of the fault surface
		# (row index = 0)
		# TODO: if point lies inside or on the surface projection of the surface boundary
		# return 0 without doing any extra calculation.
		g = Geod(ellps="sphere")
		point_list = self.fault_surf.surface.ravel()
		lons_surf = numpy.array([point_list[i].longitude for i in range(len(point_list))])
		lats_surf = numpy.array([point_list[i].latitude for i in range(len(point_list))])
		lons_point = numpy.array([point.longitude for i in range(len(point_list))])
		lats_point = numpy.array([point.latitude for i in range(len(point_list))])
		fwd_azs,back_azs,hor_dists = g.inv(lons_point,lats_point,lons_surf,lats_surf)
		return min(hor_dists * 1e-3) # to convert from m to km
		
def isOdd(i):
	return (i%2) and True or False