from math import *
from geo import *
from rup import *
import numpy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyproj

class SimpleFaultSurface:
	"""
	Class defining Simple Fault Surface.
	"""
	
	def __init__(self,fault_trace,upper_seismo_depth,lower_seismo_depth,dip,mesh_spacing = 1):
		"""
		Represents fault surface as regular (uniformly spaced) 3D mesh of points derived from:
		fault_trace: list of points representing intersection between fault surface and Earth surface
		upper_seismo_depth: minimum depth ruptures can reach (i.e. depth to fault's top edge, in km)
		lower_seismo_depth: maximum depth ruptures can reach (i.e. depth to fault's bottom edge, in km)
		dip: dip angle (in degrees)
		mesh_spacing: spacing (in km) between mesh points
		"""
		# TODO: check that all points in fault_trace have depth = 0
		# TODO: check that there are at least two points on the fault trace and they are not coincidents
		# TODO: check upper_seismo_depth >= 0
		# TODO: check dip >0 && dip <=90
		# TODO: check mesh_spacing > 0
		# TODO: check fault trace lenght is greater than fault mesh spacing
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
		
	def getStrike(self):
		"""
		Return fault strike as average azimuth of fault trace
		segments.
		"""
		return Line(self.fault_trace).getAverageAzimuth()
		
	def getDip(self):
		"""
		Return fault dip as the average dip over the fault surface mesh.
		If the surface mesh has only one location along width or one along strike
		returns dip value as passed in the constructor.
		"""
		
		if self.surface.shape[0] > 1 and self.surface.shape[1] > 1:
			return getSurfaceDip(self.surface)
		else:
			return self.dip
			
	def getSurfacePortionStrike(self,first,last_length,last_width):
		"""
		Computes representative strike value for a surface portion as defined by:
			- first: tuple (i,j) containing indexes of surface first mesh point
			- last_length: tuple (i,j) containing indexes of surface last mesh point along length
			- last_width: tuple (i,j) containing indexes of surface last mesh point along width
			If the surface portion consists of only one node along length, then the average strike
			of the horizontal segments containing the node is returned. If the fault surface mesh
			contains only one node along length, the fault strike is returned.
			If the surface portion consists of more than one node along length, then the average strike
			of the surface portion top edge is returned.
		"""
		# TODO: check that first and last_length are valid indexes
		# 0<=first[0]<=self.surface.shape[0] - 1
		# 0<=first[1]<=self.surface.shape[1] - 1
		# last_length[0]==first[0]
		# last_length[1] >= first[1]
		
		num_nodes_along_length = (last_length[1] - first[1]) + 1

		if num_nodes_along_length == 1:
			
			if self.surface.shape[1] == 1:
				return self.getStrike()
			else:
				nearest_points = self.surface[first[0],max(0, first[1] - 1):first[1] + 2]
				return Line(nearest_points.tolist()).getAverageAzimuth()
		else:
			points = self.surface[first[0],first[1]:last_length[1]+1]
			return Line(points.tolist()).getAverageAzimuth()

	def getSurfacePortionDip(self,first,last_length,last_width):
		"""
		Computes representative dip value for a surface portion as defined by:
			- first: tuple (i,j) containing indexes of surface first mesh point
			- last_length: tuple (i,j) containing indexes of surface last mesh point along length
			- last_width: tuple (i,j) containing indexes of surface last mesh point along width
		If the fault surface consists of only one node along width, the fault surface dip is returned.
		If the surface portion consists of only one node along length, then the dip is defined
		as the average dip of the mesh cells containing the uppermost node.
		If the surface portion consists of more than one node along length and width, then the dip is defined
		as the average dip of the mesh cells constituting the surface portion.
		"""
		# TODO: check that first and last_length are valid indexes
		# 0<=first[0]<=self.surface.shape[0] - 1
		# 0<=first[1]<=self.surface.shape[1] - 1
		# last_length[0]==first[0]
		# last_length[1] >= first[1]
		
		if self.surface.shape[0] == 1:
			return self.getDip()
		
		num_nodes_along_length = (last_length[1] - first[1]) + 1
		num_nodes_along_width = (last_width[0] - first[0]) + 1
		
		if num_nodes_along_length == 1:
			return getSurfaceDip(self.surface[max(0,first[0] - 1):first[0]+2,max(0, first[1] - 1):first[1] + 2])
		else:
			return getSurfaceDip(self.surface[first[0]:last_width[0]+1,first[1]:last_length[1] + 1])
			
	def getSurfacePortionCentroid(self,first,last_length,last_width):
		"""
		Computes surface portion centroid, that is its geometrical center.
		The surface portion is determined by:
		- 'first': tuple (i,j) containing indexes of surface first mesh point
		- 'last_length': tuple (i,j) containing indexes of surface last mesh point along length
		- 'last_width': tuple (i,j) containing indexes of rupture's last mesh point along width
		"""
		# if number of grid points along length and width is odd, returns
		# location corresponding to central index
		# if number of grid points along length is even and along width is
		# odd (or viceversa), returns centroid of the central segment.
		# if number of grid points along length and width is even, returns
		# centroid of the central patch
		portion = self.surface[first[0]:last_width[0]+1,first[1]:last_length[1]+1]
		num_rows = portion.shape[0]
		num_cols = portion.shape[1]
		if isOdd(num_rows) is True and isOdd(num_cols) is True:
			return portion[int((num_rows - 1) / 2),int((num_cols - 1) / 2)]
		elif isOdd(num_rows) is True and isOdd(num_cols) is not True:
			p1 = portion[int((num_rows - 1) / 2),int(num_cols / 2 - 1)]
			p2 = portion[int((num_rows - 1) / 2),int(num_cols / 2)]
			mean_lon = (p1.longitude + p2.longitude) / 2
			mean_lat = (p1.latitude + p2.latitude) / 2
			mean_depth = (p1.depth + p2.depth) / 2
			return Point(mean_lon,mean_lat,mean_depth)
		elif isOdd(num_rows) is not True and isOdd(num_cols) is True:
			p1 = portion[int(num_rows / 2 - 1),int((num_cols - 1) / 2)]
			p2 = portion[int(num_rows / 2),int((num_cols - 1)) / 2]
			mean_lon = (p1.longitude + p2.longitude) / 2
			mean_lat = (p1.latitude + p2.latitude) / 2
			mean_depth = (p1.depth + p2.depth) / 2
			return Point(mean_lon,mean_lat,mean_depth)
		else:
			p1 = portion[int(num_rows / 2 - 1),int(num_cols / 2 - 1)]
			p2 = portion[int(num_rows / 2 - 1),int(num_cols / 2)]
			p3 = portion[int(num_rows / 2),int(num_cols / 2)]
			p4 = portion[int(num_rows / 2),int(num_cols / 2 - 1)]
			mean_lon = (p1.longitude + p2.longitude + p3.longitude + p4.longitude) / 4
			mean_lat = (p1.latitude + p2.latitude + p3.latitude + p4.latitude) / 4
			mean_depth = (p1.depth + p2.depth + p3.depth + p4.depth) / 4
			return Point(mean_lon,mean_lat,mean_depth)
		
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
		self.rupture_data = self.getRuptureData()
		
	def getRuptureData(self):
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
		or as azimuth between first and last locations of fault top edge for point ruptures)
		- dip (defined as the mean dip as measured between upper and lower)
		- rake
		- tectonic region type
		- hypocenter (defined as the centroid of the rupture surface)
		- rupture surface mesh
		- rate of occurrence
		- probability of occurrence
		"""
		# extract magnitude and rate
		mag = self.rupture_data[rupt_index]['mag']
		rate = self.rupture_data[rupt_index]['rate']
		
		# extract rupture surface
		first = self.rupture_data[rupt_index]['first']
		last_length = self.rupture_data[rupt_index]['last_length']
		last_width = self.rupture_data[rupt_index]['last_width']
		rup_surf_mesh = self.fault_surf.surface[first[0]:last_width[0]+1,first[1]:last_length[1]+1]
		
		# get strike and dip
		strike = self.fault_surf.getSurfacePortionStrike(first,last_length,last_width)
		dip  = self.fault_surf.getSurfacePortionDip(first,last_length,last_width)

		# get hypocenter
		hypocenter = self.fault_surf.getSurfacePortionCentroid(first,last_length,last_width)
		
		# Poissonian probability of one or more occurrences
		probability_occurrence = 1 - exp(-rate * self.time_span)
		
		return {'magnitude':mag,'strike':strike,'dip':dip,'rake':self.rake,
				'tectonic':self.tectonic_region_type,'hypocenter':hypocenter,
				'surface':rup_surf_mesh,
				'rate':rate,'probability':probability_occurrence}
		
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
	
def getPositionVector(longitude,latitude,depth):
	"""
	Returns the position vector (in Cartesian coordinates) of
	a geographical location defined by longitude, latitude and
	depth.
	For the equation see:
	http://mathworld.wolfram.com/SphericalCoordinates.html
	Longitude and latitudes are supposed to be in degrees and 
	depth in km."""
	# TODO: check that depth is < than earth radius. Earth
	# radius should be defined in a common place.
	earth_radius = 6371 # in km
	theta = radians(longitude)
	phi = radians(90.0 - latitude)
	x = (earth_radius - depth) * cos(theta) * sin(phi)
	y = (earth_radius - depth) * sin(theta) * sin(phi)
	z = (earth_radius - depth) * cos(phi)
	return numpy.array([x,y,z])
	
def getPlaneDip(p1,p2,p3):
	"""
	Computes dip (that is inclination) of a plane with respect to
	Earth surface. The plane is defined by three corners: p1,p2,p3.
	p1 and p2 are supposed to be the upper left and right corners,
	and p3 the lower left corner.
	The dip of the plane is obtained by calculating the vector 
	normal to the plane, and the vector normal to the Earth surface
	at the upper left corner location. The angle
	between these two vectors is the dip angle.
	"""
	# for each point get corresponding cartesian coordinates
	P1 = getPositionVector(p1.longitude,p1.latitude,p1.depth)
	P2 = getPositionVector(p2.longitude,p2.latitude,p2.depth)
	P3 = getPositionVector(p3.longitude,p3.latitude,p3.depth)
	# define vectors p1p2 and p1p3, that is vector connecting upper
	# left and upper right corners and upper left and lower left corners.
	p1p2 = [P2[0]-P1[0],P2[1]-P1[1],P2[2]-P1[2]]
	p1p3 = [P3[0]-P1[0],P3[1]-P1[1],P3[2]-P1[2]]
	# compute normal vector as cross product of p1p2 and p1p3 and normalize it
	normal = numpy.cross(p1p3,p1p2)
	normal = normal / sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
	# compute unit vector normal to earth surface at p1
	normal_p1 = P1 / sqrt(P1[0]**2 + P1[1]**2 + P1[2]**2)
	# compute angle in between: this is the dip
	dip = degrees(acos(numpy.dot(normal_p1,normal)))
	return dip
	
def getSurfaceDip(surface):
	"""
	Computes dip (that is inclination) of a surface.
	The dip is computed as the average dip of the 
	mesh cells in the first row of the surface.
	The dip of each mesh cell is obtained by calculating the vector normal to each mesh
	cell, and the vector normal to the Earth surface at the cell location. The angle
	between these two vectors is the dip angle.
	The method assumes the surface has constant dip
	over depth (that's why only cells in the first
	row are considered).
	"""
	average_dip = 0.0
	for i in range(surface.shape[1] - 1):
		p1 = surface[0,i]
		p2 = surface[0,i + 1]
		p3 = surface[1,i]
		average_dip = average_dip + getPlaneDip(p1,p2,p3)
	return average_dip / (surface.shape[1] - 1)
	
def plotSimpleFaultSource(src):
	"""
	Plots simple fault source (src)
	"""

	#matplotlib.rc('axes', facecolor='k')
	# create figure
	fig = plt.figure()
	ax = fig.gca(projection='3d')

	# get orthogonal projection centered around
	# fault trace middle point
	mean_lon = (src.fault_surf.fault_trace[0].longitude + \
				src.fault_surf.fault_trace[-1].longitude) / 2
	mean_lat = (src.fault_surf.fault_trace[0].latitude + \
				src.fault_surf.fault_trace[-1].latitude) / 2
	proj = pyproj.Proj(proj='ortho', lat_0=mean_lat, lon_0=mean_lon, units='km', preserve_units=True)

	# extract fault trace coordinates
	x_top_edge,y_top_edge,depths_top_edge = getLineCartesianCoordinates(proj,Line(src.fault_surf.fault_trace))

	# extract surface points coordinates
	x_surf_points, y_surf_points, depths_surf_points = getSurfaceCartesianCoordinates(proj,src.fault_surf.surface)

	# plot fault trace
	ax.plot(x_top_edge, y_top_edge, depths_top_edge, label='Fault Trace',color='r',linewidth=3)

	# plot surface mesh points
	ax.scatter(x_surf_points, y_surf_points, depths_surf_points, c='w', marker='o')

	miny = numpy.min(y_surf_points)
	maxy = numpy.max(y_surf_points)
	minx = numpy.min(x_surf_points)
	maxx = numpy.max(x_surf_points)
	ax.set_ylim3d(miny - 1,maxy + 1)
	ax.set_xlim3d(minx - 1,maxx + 1)
	ax.set_xlabel('Along longitude (km)')
	ax.set_ylabel('Along latitude (km)')
	ax.set_zlabel('Along depth (km)')
	ax.legend()

	plt.savefig('fault_surface.png', dpi=100)
	del ax
	plt.clf()

	# then plot ruptures
	rupture_data = src.getRuptureData()
	i = 0
	for data in rupture_data:
		i += 1
		fig = plt.figure()
		ax = fig.gca(projection='3d')

		# plot fault surface points coordiantes
		ax.scatter(x_surf_points, y_surf_points, depths_surf_points, color='k', marker='.')

		# extract rupture surface points coordinates
		first = data['first']
		last_length = data['last_length']
		last_width = data['last_width']
		surf = src.fault_surf.surface[first[0]:last_width[0]+1,first[1]:last_length[1]+1]
		x_rupsurf_points, y_rupsurf_points, depths_rupsurf_points = getSurfaceCartesianCoordinates(proj,surf)

		# plot rupture surface mesh points
		ax.scatter(x_rupsurf_points, y_rupsurf_points, depths_rupsurf_points, c='m', marker='o',alpha=1.0)

		# compute rupture area
		#rup_area = getSurfacePortionArea(src.fault_surf.surface,first,last_length,last_width)
		# compute expected rupture area
		#ex_rup_area = src.mag_scaling_rel.getMedianArea(data['mag'])
		#title = 'predicted rupture area (km^2): %s\nexpected rupture area (km^2): %s' % (rup_area, ex_rup_area)
		#ax.set_title(title)

		ax.legend()
		ax.set_ylim3d(miny - 1,maxy + 1)
		ax.set_xlim3d(minx - 1,maxx + 1)
		ax.set_xlabel('Along longitude (km)')
		ax.set_ylabel('Along latitude (km)')
		ax.set_zlabel('Along depth (km)')


		filename = str('rup%s' % i) + '.png'
		plt.savefig(filename, dpi=100)
		plt.clf()
		
def getLineCartesianCoordinates(proj,line):
	lons = [p.longitude for p in line.point_list]
	lats = [p.latitude for p in line.point_list]
	depths = [-p.depth for p in line.point_list]
	x, y = proj(lons,lats)
	return x,y,depths

def getSurfaceCartesianCoordinates(proj,surf):
	"""
	surf: 2D numpy array of Locations
	"""
	surf_points = surf.flatten().tolist()
	lons_surf_points = [p.longitude for p in surf_points]
	lats_surf_points = [p.latitude for p in surf_points]
	depths_surf_points = [-p.depth for p in surf_points]
	x_surf_points, y_surf_points = proj(lons_surf_points,lats_surf_points)
	return x_surf_points, y_surf_points, depths_surf_points