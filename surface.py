import abc

class BaseSurface(object):
	"""
	Abstract class representing a 3D quadrilateral surface as 
	a 2D grid of equally spaced geographical points.
	Under the constrain of equally spaced points, the surface
	can represent a parallelogram or a rectangle.
	
	: param grid:
		2D numpy array containing grid points.
	
	: param grid_spacing:
		distance between points
	"""
	__metaclass__ = abc.ABCMeta
	
	def __init__(self,grid,grid_spacing):
		self.grid = grid
		self.grid_spacing = grid
		
	def getSurfaceCentroid(self):
		"""
		Computes surface centroid, that is its geometrical center.
		"""
		# if number of grid points along length and width is odd, returns
		# location corresponding to central index
		# if number of grid points along length is even and along width is
		# odd (or viceversa), returns centroid of the central segment.
		# if number of grid points along length and width is even, returns
		# centroid of the central patch
		num_rows = self.grid.shape[0]
		num_cols = self.grid.shape[1]
		if __isOdd(num_rows) is True and __isOdd(num_cols) is True:
			return self.grid[int((num_rows - 1) / 2),int((num_cols - 1) / 2)]
		elif __isOdd(num_rows) is True and __isOdd(num_cols) is not True:
			p1 = self.grid[int((num_rows - 1) / 2),int(num_cols / 2 - 1)]
			p2 = self.grid[int((num_rows - 1) / 2),int(num_cols / 2)]
			mean_lon = (p1.longitude + p2.longitude) / 2
			mean_lat = (p1.latitude + p2.latitude) / 2
			mean_depth = (p1.depth + p2.depth) / 2
			return Point(mean_lon,mean_lat,mean_depth)
		elif __isOdd(num_rows) is not True and __isOdd(num_cols) is True:
			p1 = self.grid[int(num_rows / 2 - 1),int((num_cols - 1) / 2)]
			p2 = self.grid[int(num_rows / 2),int((num_cols - 1)) / 2]
			mean_lon = (p1.longitude + p2.longitude) / 2
			mean_lat = (p1.latitude + p2.latitude) / 2
			mean_depth = (p1.depth + p2.depth) / 2
			return Point(mean_lon,mean_lat,mean_depth)
		else:
			p1 = self.grid[int(num_rows / 2 - 1),int(num_cols / 2 - 1)]
			p2 = self.grid[int(num_rows / 2 - 1),int(num_cols / 2)]
			p3 = self.grid[int(num_rows / 2),int(num_cols / 2)]
			p4 = self.grid[int(num_rows / 2),int(num_cols / 2 - 1)]
			mean_lon = (p1.longitude + p2.longitude + p3.longitude + p4.longitude) / 4
			mean_lat = (p1.latitude + p2.latitude + p3.latitude + p4.latitude) / 4
			mean_depth = (p1.depth + p2.depth + p3.depth + p4.depth) / 4
			return Point(mean_lon,mean_lat,mean_depth)

	def getShortestDistance(self,point):
		"""
		Computes shortest distance (in km) from surface to point. The shortest distance 
		is computed as the minimum among the distances between the point of 
		interest and the points constituting the grid.
		"""
		g = Geod(ellps="sphere")
		point_list = self.grid.ravel()
		lons_surf = numpy.array([point_list[i].longitude for i in range(len(point_list))])
		lats_surf = numpy.array([point_list[i].latitude for i in range(len(point_list))])
		lons_point = numpy.array([point.longitude for i in range(len(point_list))])
		lats_point = numpy.array([point.latitude for i in range(len(point_list))])
		fwd_azs,back_azs,hor_dists = g.inv(lons_point,lats_point,lons_surf,lats_surf)
		hor_dists = hor_dists * 1e-3 # to convert from m to km
		vert_dists = numpy.array([(point.depth - point_list[i].depth) for i in range(len(point_list))])
		dists = numpy.sqrt(hor_dists**2 + vert_dists**2)
		return min(dists)
		
	def getShortestDistanceFromSurfaceProjection(self,point):
		"""
		Computes shortest distance from Earth's surface projection of the surface
		to point. The shortest distance is computed as the minimum among the 
		distances between the point of interest and the surface projections 
		of the points constituting the grid.
		"""
		# TODO: if point lies inside or on the surface projection of the surface boundary
		# return 0 without doing any extra calculation.
		# TODO: if rupture is vertical compute distance only to points on the top edge of 
		# the surface
		g = Geod(ellps="sphere")
		point_list = self.grid.ravel()
		lons_surf = numpy.array([point_list[i].longitude for i in range(len(point_list))])
		lats_surf = numpy.array([point_list[i].latitude for i in range(len(point_list))])
		lons_point = numpy.array([point.longitude for i in range(len(point_list))])
		lats_point = numpy.array([point.latitude for i in range(len(point_list))])
		fwd_azs,back_azs,hor_dists = g.inv(lons_point,lats_point,lons_surf,lats_surf)
		return min(hor_dists * 1e-3) # to convert from m to km
		
	def getDistanceX(self,point):
		"""
		Computes horizontal distance (in km) to the surface top edge
		measured perpendicular to the surface azimuth (computed
		considering first and last point in the surface top edge.)
		"""
		# TODO: extends surface top edge, along the azimuth direction
		
	def __isOdd(i):
		return (i%2) and True or False