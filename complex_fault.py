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
		- the firzt point in the fault intermediate edge lies on the left side of the fault surface
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
	
	def __init__(self,fault_surf,freq_mag_dist,mag_scaling_rel,rake,rup_aspect_ratio,tectonic_region_type,time_span,area_tol):
		"""
		fault_surf: complex fault surface
		freq_mag_dist: frequency magnitude distribution
		mag_scaling_rel: magnitude scaling relationship
		rake: rake angle (-180 <= rake <= 180)
		rup_aspect_ratio: rupture aspect ratio (> 0)
		tectonic_region_type: tectonic region type
		time_span: time span (>= 0)
		tol = tolerance level (in percent)
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
		self.area_tol = area_tol
		
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
		
		# TODO: maybe self.area_tol = 100 * numpy.max(cells_area) / ex_rup_area
		
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
				
			else:
			
				# loop over ruptures' upper left corners
				for i in range(self.fault_surf.surface.shape[0] - 1):
					for j in range(self.fault_surf.surface.shape[1] - 1):
					
						# compute possible rupture areas, starting from node
						# (i,j) by accumulating cells areas along length (i.e. rows)
						# and along width (i.e. columns)
						rup_areas = numpy.add.accumulate(cells_area[i:,j:],axis=1)
						rup_areas = numpy.add.accumulate(rup_areas,axis=0)
						
						# compute corresponding rupture aspect ratios,
						# starting from node (i,j)
						rup_lengths = numpy.add.accumulate(cells_lengths[i:,j:],axis=1)
						rup_widths = numpy.add.accumulate(cells_widths[i:,j:],axis=0)
						aspect_ratios = rup_lengths / rup_widths
						
						# extract node indexes giving rupture areas
						# close to expected rupture area (within tolerance)
						rup_indexes_area = numpy.where(100 * abs(rup_areas - ex_rup_area) / ex_rup_area<= self.area_tol)
						rup_indexes_area_z = zip(rup_indexes_area[0],rup_indexes_area[1])
						
						# if there are not just continue (that is, in the current position
						# the fault surface cannot accomodate a rupture with the
						# expected area [within tolerance])
						if len(rup_indexes_area[0]) == 0:
							continue
						
						# among the ruptures consistent with the expected rupture area,
						# extract the one that has the aspect ratio closest to the one
						# given in the constructor
						rup_index = numpy.where(abs(aspect_ratios - self.rup_aspect_ratio) == numpy.min(abs(aspect_ratios[rup_indexes_area] - self.rup_aspect_ratio)))
						
						# extract last nodes along length and width
						# the plus 1 is due to the fact that rup_index
						# corresponds to the cell index, while
						# we are interested in the surface last node index
						last_length = (i,j + rup_index[1][0] + 1)
						last_width = (i + rup_index[0][0] + 1,j)
								
						data = {'mag':mag,
								'rate':rate,
								'first':(i,j),
								'last_length':last_length,
								'last_width':last_width}
						rupture_data.append(data)
						
						# if the rupture touches the right boundary of the fault, break,
						# that is continue on the next row
						if last_length[1] == self.fault_surf.surface.shape[1] - 1:
							break
						
		return rupture_data
		
		
		
		#fault_surface_area = getSurfacePortionArea(self.fault_surf.surface,(0,0),(0,self.fault_surf.surface.shape[1] - 1),(self.fault_surf.surface.shape[0]-1,0))
		
		# loop over magnitude values, For each magnitude compute expected rupture area.
		# If expected rupture area is greater than fault surface area, return fault surface
		# as rupture surface.
		# If expected rupture area is smaller than fault surface area, loop over fault source
		# nodes (except nodes on the fault bottom and right edges).
		# For each node (for the purpose of the algorithm identified by the name nucleation node),
		# determine rupture surface last node along length and along width.
		# Along length and width last nodes are determined by finding the nodes
		# (along rupture horizontal and vertical edges) that defines a rupture surface with an area
		# that is closest to the expected rupture area.
		# if the rupture surface node opposite to the nucleation node, lies on a fault surface edge
		# then the rupture is considered only if the percentage difference between rupture surface
		# area and expted rupture area is smaller than a tolerance level.
		#for mag,rate in occurrence_rates:
			
			# compute expected rupture surface area
			#ex_rup_area = self.mag_scaling_rel.getMedianArea(mag)
			#ex_rup_length = sqrt(ex_rup_area * self.rup_aspect_ratio)
			#ex_rup_width = ex_rup_area / ex_rup_length		
			
			#if ex_rup_area >= fault_surface_area:
				# return indexes corresponding to the entire surface
			#	data = {'mag':mag,
			#			'rate':rate,
			#			'first':(0,0),
			#			'last_length':(0,self.fault_surf.surface.shape[1] - 1),
			#			'last_width':(self.fault_surf.surface.shape[0]-1,0)}
			#	rupture_data.append(data)
			#else:
			#	# loop over ruptures' upper left corners
			#	for i in range(self.fault_surf.surface.shape[0] - 1):
			#		for j in range(self.fault_surf.surface.shape[1] - 1):
			#			
			#			# extract lines corresponding to length and width
			#			line_l = Line(self.fault_surf.surface[i,j:].tolist())
			#			line_w = Line(self.fault_surf.surface[i:,j].tolist())
						
						# compute length from current upper left corner of
						# the rupture 
			#			max_length = line_l.getLength()
						
						# compute rupture length corresponding to rupture area closest to expected rupture area
			#			print 'optimizing...'
			#			args = (self.fault_surf.surface,(i,j),self.rup_aspect_ratio,ex_rup_area)			
			#			xtol = self.fault_surf.mesh_spacing/2
			#			maxfun = len(line_l.point_list)/2
			#			length= optimize.fminbound(ruptureAreaRelativeDiff, 0, max_length, args=args, xtol=xtol, maxfun=5, full_output=0, disp=0)
						
						# get points corresponding to length and width as obtained from length / aspect_ratio
			#			last_length = numpy.where(self.fault_surf.surface == line_l.getClosestPoint(length))
			#			last_width = numpy.where(self.fault_surf.surface == line_w.getClosestPoint(length/self.rup_aspect_ratio))
			#			last_length  = (last_length[0][0],last_length[1][0])
			#			last_width = (last_width[0][0],last_width[1][0])
						
			#			print 'computing area...'
						# compute rupture surface area
			#			comp_rup_area = getSurfacePortionArea(self.fault_surf.surface,(i,j),last_length,last_width)
						
						# compute percentage difference between expected and computed area
						# if higher than tolerance, break
			#			percent_relative_diff = 100 * abs(comp_rup_area - ex_rup_area) / ex_rup_area
			#			print (i,j)
			#			print 'expected rupture area: %s, computed rupture area: %s, percent diff: %s' % (ex_rup_area,comp_rup_area,percent_relative_diff)
			#			if last_width[0] == self.fault_surf.surface.shape[0] - 1 or last_length[1] == self.fault_surf.surface.shape[1] - 1:
			#				if percent_relative_diff > self.tol:
			#					break

def ruptureAreaRelativeDiff(length,surface,first,aspect_ratio,expected_area):
	"""
	Computes relative percent difference between 'expected_area', and
	surface portion area defined over 'surface', and identified
	by:
	- length: portion length (along top edge)
	- first: (i,j) tuple identifing portion upper left corner node
	- aspect_ratio: (top edge) length / (left edge) width
	"""
	
	# computes surface portion width based on length and aspect ratio
	width = length / aspect_ratio
	
	# extract lines corresponding to length and width
	line_l = Line(surface[first[0],first[1]:].tolist())
	line_w = Line(surface[first[0]:,first[1]].tolist())
	
	# get closest point along length and with
	# and compute corresponding area
	last_length = numpy.where(surface == line_l.getClosestPoint(length))
	last_width = numpy.where(surface == line_w.getClosestPoint(width))
	last_length  = (last_length[0][0],last_length[1][0])
	last_width = (last_width[0][0],last_width[1][0])
	computed_area =  getSurfacePortionArea(surface,first,last_length,last_width)
	
	# return (percent) relative differences
	f =  100 * abs(computed_area - expected_area) / expected_area
	
	return f
						
def getSurfacePortionArea(surface,first,last_length,last_width):
	"""
	Computes surface portion area.
	The surface portion is determined by:
	- 'first': tuple (i,j) containing indexes of surface first mesh point
	- 'last_length': tuple (i,j) containing indexes of surface last mesh point along length
	- 'last_width': tuple (i,j) containing indexes of rupture's last mesh point along width
	"""

	# TODO: check that first,last_length,last_width are valid indexes
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

		ax.set_xlabel('Along longitude (km)')
		ax.set_ylabel('Along latitude (km)')
		ax.set_zlabel('Along depth (km)')
		ax.legend()
		
		filename = str('rup%s' % i) + '.png'
		plt.savefig(filename, dpi=100)
		del ax
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