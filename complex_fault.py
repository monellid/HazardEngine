from math import *
from geo import *
from rup import *
import numpy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyproj

class ComplexFaultSurface:
	"""
	Class defining Complex Fault Surface.
	"""
	def __init__(self,fault_top_edge,fault_bottom_edge,mesh_spacing,fault_intermediate_edge = None):
		"""
		Represents fault surface as 3D mesh of points obtained from:
		fault_top_edge: list of points representing fault top edge
		fault_intermediate_edge: list of points representing fault intermediate edge (default None)
		fault_bottom_edge: list of points representing fault bottom edge
		mesh_spacing: average spacing between grid nodes
		The class assumes that with respect to the surface centroid:
		- the first point in the fault top edge is the fault upper left corner
		- the last point in the fault top edge is the fault upper right corner
		- the first point in the fault bottom edge is the fault lower left corner
		- the last point in the fault bottom edge is the fault lower right corner
		- the first point in the fault intermediate edge lies on the left side of the fault surface
		- the last point in the fault intermediate edge lies on the right side of the fault surface
		"""
		#TODO: check fault_top_edge != fault_bottom_edge != fault_intermediate_edge (i.e. they cannot contain the same set of points)
		#TODO: check fault_top_edge, fault_bottom_edge, and fault_intermediate_edge do not share any points
		#TODO: check that by joining fault corners the resulting
		# polygon is valid. That is defines two lines connecting upper and lower left corner, and
		# upper and lower right corner. Check that the two lines are not intersecting each other
		# on the lat, lon plane.
		# do the same between fault_top_edge and fault_intermediate_edge and fault_intermediate_edge and fault_bottom_edge
		#TODO: check mesh_spacing > 0
		self.fault_top_edge = fault_top_edge
		self.fault_bottom_edge = fault_bottom_edge
		self.mesh_spacing = mesh_spacing
		self.fault_intermediate_edge = fault_intermediate_edge
		self.surface = self.getSurfaceMesh()
		
	def getSurfaceMesh(self):
		"""
		Computes 3D mesh representing fault surface.
		"""
		# computes mean fault edge length
		length_top_edge = self.fault_top_edge.getLength()
		length_bottom_edge = self.fault_bottom_edge.getLength()
		if self.fault_intermediate_edge is not None:
			length_intermediate_edge = self.fault_intermediate_edge.getLength()
			mean_length = (length_top_edge + length_bottom_edge + length_intermediate_edge) / 3
		else:
			mean_length = (length_top_edge + length_bottom_edge) / 2
			
		# compute number of points on each edge based on mean length and mesh_spacing
		num_points_length = int(round(mean_length / self.mesh_spacing)) + 1
		
		# resample edges in num_nodes points
		top_edge = self.fault_top_edge.getResampledLineInNpoints(num_points_length)
		bottom_edge = self.fault_bottom_edge.getResampledLineInNpoints(num_points_length)
		if self.fault_intermediate_edge is not None:
			intermediate_edge = self.fault_intermediate_edge.getResampledLineInNpoints(num_points_length)
			
		# compute mean fault width
		average_width = 0.0
		for i, point in enumerate(top_edge.point_list):
			
			point_list = []
			point_list.append(point)
			if self.fault_intermediate_edge is not None:
				point_list.append(intermediate_edge.point_list[i])
			point_list.append(bottom_edge.point_list[i])
			
			average_width = average_width + Line(point_list).getLength()
			
		average_width = average_width / num_points_length
		
		# compute number of points along width, based on average width and mesh_spacing
		num_points_width = int(round(average_width / self.mesh_spacing)) + 1 
			
		# create surface
		# loop over nodes in the upper edge, define line connecting upper node with
		# bottom node (passing through intermediate node if defined)
		# then resample line in num_points
		mesh_points = []
		for i, point in enumerate(top_edge.point_list):
			
			point_list = []
			point_list.append(point)
			if self.fault_intermediate_edge is not None:
				point_list.append(intermediate_edge.point_list[i])
			point_list.append(bottom_edge.point_list[i])
			points = Line(point_list).getResampledLineInNpoints(num_points_width)
			mesh_points.extend(points.point_list)
			
		# organize mesh points into a 2D array
		# number of rows corresponds to number of points along dip
		# number of columns corresponds to number of points along strike
		surface = numpy.array(mesh_points)
		surface = surface.reshape(num_points_length,num_points_width)
		surface = numpy.transpose(surface)
		return surface
		
class ComplexFaultSource:
	"""
	Class defining complex fault source.
	"""
	
	def __init__(self,fault_surf,freq_mag_dist,mag_scaling_rel,rake,rup_aspect_ratio,tectonic_region_type,time_span):
		"""
		fault_surf: complex fault surface
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
		NOTE: the algorithm is designed in a way that each rupture is defined by at least two points
		along length, and two points along width (independently of magnitude and mesh spacing).
		That is mesh_spacing * mesh_spacing is the 'area quantum'. All ruptures with area smaller than
		the 'area quantum' are modelled as having an area equal to the 'area quantum'. 
		"""
		
		rupture_data = []
		
		# compute fault surface area
		fault_surface_area = getSurfacePortionArea(self.fault_surf.surface,(0,0),(0,self.fault_surf.surface.shape[1] - 1),(self.fault_surf.surface.shape[0]-1,0))
		
		# computes surface mesh cells areas, lengths and widths
		cells_area = numpy.ndarray([self.fault_surf.surface.shape[0]-1,self.fault_surf.surface.shape[1]-1])
		cells_lengths = numpy.ndarray([self.fault_surf.surface.shape[0]-1,self.fault_surf.surface.shape[1]-1])
		cells_widths = numpy.ndarray([self.fault_surf.surface.shape[0]-1,self.fault_surf.surface.shape[1]-1])
		for i in range(self.fault_surf.surface.shape[0]-1):
			for j in range(self.fault_surf.surface.shape[1]-1):
				cells_area[i,j] = getSurfacePortionArea(self.fault_surf.surface,(i,j),(i,j+1),(i+1,j))
				cells_lengths[i,j] = Line(self.fault_surf.surface[i,j:j+2].tolist()).getLength()
				cells_widths[i,j] = Line(self.fault_surf.surface[i:i+2,j].tolist()).getLength()
		
		occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
		
		# normalization factors. Dictionary containing, per each magnitude value,
		# the number of ruptures defined.
		norm_f = {}
		
		for mag,rate in occurrence_rates:
			
			# compute expected rupture surface area, length and width
			ex_rup_area = self.mag_scaling_rel.getMedianArea(mag)
			ex_rup_length = sqrt(ex_rup_area * self.rup_aspect_ratio)
			ex_rup_width = ex_rup_area / ex_rup_length
			
			if ex_rup_area >= fault_surface_area:
				# return indexes corresponding to the entire surface
				data = {'mag':mag,
						'rate':rate,
						'first':(0,0),
						'last_length':(0,self.fault_surf.surface.shape[1] - 1),
						'last_width':(self.fault_surf.surface.shape[0]-1,0)}
				rupture_data.append(data)
				norm_f[mag] = 1
			else:
			
				count  = 0
				# array containing boolean values (set to False), used to check
				# if a node at the fault bottom edge has been already used as
				# bottom left corner of a rupture.
				visited_nodes = numpy.ndarray((self.fault_surf.surface.shape[1]),dtype=bool)
				visited_nodes[:] = False
				# loop over ruptures' upper left corners
				for i in range(self.fault_surf.surface.shape[0] - 1):
					for j in range(self.fault_surf.surface.shape[1] - 1):
						
						# compute possible rupture lengths from current node
						rup_lengths = numpy.add.accumulate(cells_lengths[i:,j:],axis=1)
						
						# extract the node that corresponds to a length closest to
						# the expected one
						last_length_idx = numpy.where(abs(rup_lengths[0,:] - ex_rup_length) == numpy.min(abs(rup_lengths[0,:] - ex_rup_length)))
						
						# compute possible rupture areas, starting from node
						# (i,j) by accumulating cells areas along length (i.e. rows)
						# and along width (i.e. columns)
						rup_areas = numpy.add.accumulate(cells_area[i:,j:],axis=1)
						rup_areas = numpy.add.accumulate(rup_areas,axis=0)
						
						# extract node corresponding to a width giving area closest to the expected one
						# but keeping the along length node fixed
						last_width_idx = numpy.where(abs(rup_areas[:,last_length_idx[0][0]] - ex_rup_area) == numpy.min(abs(rup_areas[:,last_length_idx[0][0]] - ex_rup_area)))
						
						# if the last node along width lies on the fault bottom boundary
						# adjust rupture length so as to reach the optimal area value
						# that is, sacrifice rupture aspect ratio for rupture area
						if i + last_width_idx[0][0] + 1 == self.fault_surf.surface.shape[0] - 1:
							last_length_idx = numpy.where(abs(rup_areas[last_width_idx[0][0],:] - ex_rup_area) == numpy.min(abs(rup_areas[last_width_idx[0][0],:] - ex_rup_area)))
						
						# if the last node along width lies on the fault bottom boundary, check if
						# it has been already used as "bottom left corner" for a previous rupture.
						# If not consider the rupture (and set visited_node to True), otherwise continue.
						if i + last_width_idx[0][0] + 1 == self.fault_surf.surface.shape[0] - 1:
							if visited_nodes[j] == False:
								visited_nodes[j] = True
							else:
								continue
						
						# extract last nodes along length and width
						# the plus 1 is due to the fact that rup_index
						# corresponds to the cell index, while
						# we are interested in the surface last node index
						last_length = (i,j + last_length_idx[0][0] + 1)
						last_width = (i + last_width_idx[0][0] + 1,j)
								
						data = {'mag':mag,
								'rate':rate,
								'first':(i,j),
								'last_length':last_length,
								'last_width':last_width}
						rupture_data.append(data)
						count = count + 1
						
						# if the rupture touches the right boundary of the fault, break,
						# that is continue on the next row
						if last_length[1] == self.fault_surf.surface.shape[1] - 1:
							break
							
					# we stop at the first rupture touching both the bottom and right fault
					# boundaries
					if rupture_data[-1]['last_length'][1] == self.fault_surf.surface.shape[1] - 1 and rupture_data[-1]['last_width'][0] == self.fault_surf.surface.shape[0] - 1:
						break
				
				# store how many ruptures have been defined for the given magnitude	
				norm_f[mag] = count
						
		# loop over defined ruptures and scale rates
		for i in range(len(rupture_data)):
			rupture_data[i]['rate'] = rupture_data[i]['rate'] / norm_f[rupture_data[i]['mag']]
			
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
		- strike
		- dip
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
		strike = getSurfacePortionStrike(self.fault_surf.surface,first,last_length,last_width)
		dip  = getSurfacePortionDip(self.fault_surf.surface,first,last_length,last_width)

		# get hypocenter
		hypocenter = getSurfacePortionCentroid(self.fault_surf.surface,first,last_length,last_width)

		# Poissonian probability of one or more occurrences
		probability_occurrence = 1 - exp(-rate * self.time_span)

		return {'magnitude':mag,'strike':strike,'dip':dip,'rake':self.rake,
				'tectonic':self.tectonic_region_type,'hypocenter':hypocenter,
				'surface':rup_surf_mesh,
				'rate':rate,'probability':probability_occurrence}
				
def getSurfacePortionCentroid(surface,first,last_length,last_width):
	"""
	Computes surface portion centroid. 
	The algorithm works as follows. The surface portion is splitted into
	triangular facets. The centroid is then computed as the center of mass
	of the surgface portion,following the equation in 
	http://mathworld.wolfram.com/GeometricCentroid.html
	where the 'mass' of each triangle is its area.
	Return centroid as numpy.array([longitude,latitude,depth])
	"""
	x = []
	y = []
	z = []
	area = []
	for i in range(first[0],last_width[0]):
		for j in range(first[1],last_length[1]):

			# define the two triangles constituting the mesh cell
			# lower triangle
			t1 = Triangle(surface[i,j],surface[i+1,j],surface[i+1,j+1])
			# upper triangle
			t2 = Triangle(surface[i,j],surface[i,j+1],surface[i+1,j+1])
			
			area1 = t1.getArea()
			centroid1 = t1.getCentroid()
			
			area2 = t2.getArea()
			centroid2 = t2.getCentroid()
			
			x.append(centroid1[0])
			x.append(centroid2[0])
			
			y.append(centroid1[1])
			y.append(centroid2[1])
			
			z.append(centroid1[2])
			z.append(centroid2[2])
			
			area.append(area1)
			area.append(area2)
	
	x_centroid = numpy.sum(numpy.array(area) * numpy.array(x)) / numpy.sum(numpy.area)
	y_centroid = numpy.sum(numpy.array(area) * numpy.array(y)) / numpy.sum(numpy.area)
	z_centroid = numpy.sum(numpy.array(area) * numpy.array(z)) / numpy.sum(numpy.area)
	
	# convert to spherical coordinates
	return getSphericalPositionVector(x_centroid,y_centroid,z_centroid)
		
def getSurfacePortionStrike(surface,first,last_length,last_width):
	"""
	Computes surface portion strike. Assumes the surface to be defined by a 2D numpy array
	of locations, with variable spacing. Locations on different rows are not requested to
	be aligned along the same direction.
	The algorithm computes the surface portion strike, by calculating first the average azimuth
	of each surface line (that is locations defined in the same row), and then the weigthed average azimuth
	of the surface lines' average azimuths (the weights are the surface lines lengths).
	NOTE: given that this method is used for the complex fault source, each surface portion
	is defined by at least two points along lenght and along width, so there is no need
	to deal with the case of a surface portion defined by a single point along length and/or
	along width.
	"""
	azimuths = []
	weights = []
	for i in range(first[0],last_width[0] + 1):
		line = Line(surface[i,first[1]:last_length[1] + 1])
		azimuths.append(radians(line.getAverageAzimuth()))
		weights.append(line.getLength())
		
	# convert from polar to cartesian coordinates
	vectors = []
	for i in range(len(azimuths)):
		vectors.append(numpy.array([weights[i] * sin(azimuths[i]),weights[i] * cos(azimuths[i])]))
		
	# sum all vectors. this represents the mean direction,
	# from which we can extract the mean angle
	v = vectors[0]
	for i in range(1,len(vectors)):
		v = v + vectors[i]
		
	# extract angle
	strike = degrees(atan(v[0] / v[1]))
	
	if strike < 0:
		strike = strike + 360.0
		
	if strike >= 360.0:
		strike = strike - 360.0
	
	return strike
	
def getSurfacePortionDip(surface,first,last_length,last_width):
	"""
	Computes surface portion dip. Assumes the surface to be defined by a 2D numpy array
	of locations, with variable spacing. Locations on different rows are not requested to
	be aligned along the same direction.
	The algorithm computes the surface portion dip by calculating the weighted average of the
	surface portion mesh cells (the weights being the mesh cells). The dip of each mesh 
	cell is calculated by splitting the cell into two triangles, and calculating the 
	dip of each triangle (each triangle defines a plane with a certain inclination/dip 
	with respect to the Earth surface). The dip of each cell is the weighted average 
	of the triangles' dips (weights are the triangles' areas).
	"""
	cell_dips = []
	cell_areas = []
	for i in range(first[0],last_width[0]):
		for j in range(first[1],last_length[1]):

			# define the two triangles constituting the mesh cell
			# lower triangle
			t1 = Triangle(surface[i,j],surface[i+1,j],surface[i+1,j+1])
			# upper triangle
			t2 = Triangle(surface[i,j],surface[i,j+1],surface[i+1,j+1])
			
			# compute dips and areas of the two triangles
			dip1 = t1.getInclination()
			area1 = t1.getArea()
			dip2 = t2.getInclination()
			area2 = t2.getArea()
			
			cell_dips.append((area1 * dip1 + area2 * dip2) / (area1 + area2))
			cell_areas.append(area1 + area2)
			
	# compute weighted mean
	return numpy.sum(numpy.array(cell_dips) * numpy.array(cell_areas)) / numpy.sum(numpy.array(cell_areas))
	
def getSurfacePortionArea(surface,first,last_length,last_width):
	"""
	Computes surface portion area. Assumes the surface to be defined by a 2D numpy array
	of locations, with variable spacing. Locations on different rows are not requested to
	be aligned along the same direction.
	The algorihtm computes the surface portion area, by summing up the areas of the mesh
	cells composing the surface portion. The area of each mesh cell is calculated by splitting
	the mesh cell into two triangles, and calculating the area of each triangle. The sum of the
	triangles' areas gives the mesh cell area.
	The surface portion is determined by:
	- 'first': tuple (i,j) containing indexes of surface first mesh point
	- 'last_length': tuple (i,j) containing indexes of surface last mesh point along length
	- 'last_width': tuple (i,j) containing indexes of rupture's last mesh point along width
	"""
	area = 0.0
	for i in range(first[0],last_width[0]):
		for j in range(first[1],last_length[1]):

			# define the two triangles constituting the mesh cell
			# lower triangle
			t1 = Triangle(surface[i,j],surface[i+1,j],surface[i+1,j+1])
			# upper triangle
			t2 = Triangle(surface[i,j],surface[i,j+1],surface[i+1,j+1])

			# compute area of the two triangles and sum them up
			area_t1 = t1.getArea()
			area_t2 = t2.getArea()
			area = area + area_t1 + area_t2

	return area
	
def plotComplexFaultSource(src):
	"""
	Plots complex fault source (src)
	"""
	
	#matplotlib.rc('axes', facecolor='k')
	# create figure
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	
	# get orthogonal projection centered around
	# fault top trace middle point
	mean_lon = (src.fault_surf.fault_top_edge.point_list[0].longitude + \
				src.fault_surf.fault_top_edge.point_list[-1].longitude) / 2
	mean_lat = (src.fault_surf.fault_top_edge.point_list[0].latitude + \
				src.fault_surf.fault_top_edge.point_list[-1].latitude) / 2
	proj = pyproj.Proj(proj='ortho', lat_0=mean_lat, lon_0=mean_lon, units='km', preserve_units=True)
	
	# extract fault and bottom edges coordinates
	x_top_edge,y_top_edge,depths_top_edge = getLineCartesianCoordinates(proj,src.fault_surf.fault_top_edge)
	x_bottom_edge, y_bottom_edge, depths_bottom_edge = getLineCartesianCoordinates(proj,src.fault_surf.fault_bottom_edge)
	
	# extract fault intermediate edge coordinates (if defined)
	if src.fault_surf.fault_intermediate_edge is not None:
		x_intermediate_edge, y_intermediate_edge, depths_intermediate_edge = \
			getLineCartesianCoordinates(proj,src.fault_surf.fault_intermediate_edge)
		
	# extract surface points coordinates
	x_surf_points, y_surf_points, depths_surf_points = getSurfaceCartesianCoordinates(proj,src.fault_surf.surface)
	
	# plot top edge
	ax.plot(x_top_edge, y_top_edge, depths_top_edge, label='Top edge',color='r',linewidth=3)
	# plot intermediate edge
	if src.fault_surf.fault_intermediate_edge is not None:
		ax.plot(x_intermediate_edge, y_intermediate_edge, depths_intermediate_edge, label='Intermediate edge',color='g',linewidth=3)
	# plot bottom edge
	ax.plot(x_bottom_edge, y_bottom_edge, depths_bottom_edge, label='Bottom edge',color='b',linewidth=3)
	
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
		rup_area = getSurfacePortionArea(src.fault_surf.surface,first,last_length,last_width)
		# compute expected rupture area
		ex_rup_area = src.mag_scaling_rel.getMedianArea(data['mag'])
		title = 'predicted rupture area (km^2): %s\nexpected rupture area (km^2): %s' % (rup_area, ex_rup_area)
		ax.set_title(title)
		
		#ax.set_aspect('equal')
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