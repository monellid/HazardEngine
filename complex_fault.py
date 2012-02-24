from math import *
from geo import *
from rup import *
import numpy
from scipy import optimize


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
	
	def __init__(self,fault_surf,freq_mag_dist,mag_scaling_rel,rake,rup_aspect_ratio,tectonic_region_type,time_span,tol):
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
		self.tol = tol
		
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
		
		fault_surface_area = getSurfacePortionArea(self.fault_surf.surface,(0,0),(0,self.fault_surf.surface.shape[1] - 1),(self.fault_surf.surface.shape[0]-1,0))
			
		occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
		
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
		for mag,rate in occurrence_rates:
			
			# compute expected rupture surface area
			ex_rup_area = self.mag_scaling_rel.getMedianArea(mag)
			
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
						
						# extract lines corresponding to length and width
						line_l = Line(self.fault_surf.surface[i,j:].tolist())
						line_w = Line(self.fault_surf.surface[i:,j].tolist())
						
						# compute length from current upper left corner of
						# the rupture 
						max_length = line_l.getLength()
						
						# compute rupture length corresponding to rupture area closest to expected rupture area
						print 'optimizing...'
						args = (self.fault_surf.surface,(i,j),self.rup_aspect_ratio,ex_rup_area)			
						xtol = self.fault_surf.mesh_spacing/2
						length, fval, ierr, numfunc = optimize.fminbound(ruptureAreaRelativeDiff, 0.0, max_length, args=args, xtol=xtol, maxfun=500, full_output=1, disp=0)
						
						# get points corresponding to length and width as obtained from length / aspect_ratio
						last_length = numpy.where(self.fault_surf.surface == line_l.getClosestPoint(length))
						last_width = numpy.where(self.fault_surf.surface == line_w.getClosestPoint(length/self.rup_aspect_ratio))
						last_length  = (last_length[0][0],last_length[1][0])
						last_width = (last_width[0][0],last_width[1][0])
						
						print 'computing area...'
						# compute rupture surface area
						comp_rup_area = getSurfacePortionArea(self.fault_surf.surface,(i,j),last_length,last_width)
						
						# compute percentage difference between expected and computed area
						# if higher than tolerance, break
						percent_relative_diff = 100 * abs(comp_rup_area - ex_rup_area) / ex_rup_area
						print (i,j)
						print 'expected rupture area: %s, computed rupture area: %s, percent diff: %s' % (ex_rup_area,comp_rup_area,percent_relative_diff)
						if last_width[0] == self.fault_surf.surface.shape[0] - 1 or last_length[1] == self.fault_surf.surface.shape[1] - 1:
							if percent_relative_diff > self.tol:
								break

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
	return 100 * abs(computed_area - expected_area) / expected_area
						
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