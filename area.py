from geo import *
from pyproj import Geod
from shapely.geometry import Polygon
from shapely.geometry import Point as shapelyPoint
from math import *
from point import *

# TODO: this should be defined in a common place.
# all distance calculations in the engine
# must use the same Geod object based on the same ellipsoid
g = Geod(ellps='sphere')
# spacing used to resample area boundary (in km).
# this value corresponds to about 1 degree at the equator
GREAT_CIRCLE_DELTA = 100.0

class PoissonianAreaSource:
	"""
	Class defining Poissonian Area source. Earthquake ruptures are distributed
	uniformly over a geographical region defined by a polygonal area. Ruptures have
	poissonian probability of occurrence.
	Given a polygon defining a geographical region on the Earth surface, the polygon
	area is discretized with a uniform grid. A point source is then defined at each
	grid point. The magnitude frequency distribution associated to each point source
	is obtained by scaling the original magnitude frequency distribution associated
	to the area by a factor equal to the area of the grid cell divided by the total
	area of the polygon.
	"""
	
	def __init__(self,area_boundary,freq_mag_dist,nodal_plane_pmf,hypo_depth_pmf,upper_seismo_depth,lower_seismo_depth,mag_scaling_rel,rupt_aspect_ratio,tectonic_region_type,time_span,grid_spacing,mesh_spacing):
		"""
		area_boundary: list of points defining a polygon on the earth surface [(lon1,lat1),(lon2,lat2),...,(lonn,latn)]
		freq_mag_dist: frequency magnitude distribution
		nodal_plane_dist: nodal plane (strike,dip,rake) probability mass function [(strike1,dip1,rake1,weight1),...,(striken,dipn,raken,weightn)]
		hypo_depth_pmf: hypocentral depth (in km) probability mass function [(hypo1,weight1),...,(hypo2,weight2)]
		upper_seismo_depth: upper seismogenic depth (in km)
		lower_seismo_depth: lower seismogenic depth (in km)
		mag_scaling_rel: magnitude-area scaling relationship
		rupt_aspect_ratio: rupture aspect ratio
		tectonic_region_type: tectonic region type
		time_span: time span (years) for calculating probability of occurrence
		grid_spacing: area grid spacing (in degrees)
		mesh_spacing: rupture surface mesh spacing (in km)
		"""
		#TODO: check polygon is well defined (no intermediate duplicated points, poligon boundary does not intersect itself)
		#TODO: in the polygon the last point must match the first point. If it's not add
		# a copy of the first point at the end of the list.
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
		self.area_boundary = area_boundary
		self.freq_mag_dist = freq_mag_dist
		self.nodal_plane_pmf = nodal_plane_pmf
		self.hypo_depth_pmf = hypo_depth_pmf
		self.upper_seismo_depth = upper_seismo_depth
		self.lower_seismo_depth = lower_seismo_depth
		self.mag_scaling_rel = mag_scaling_rel
		self.rupt_aspect_ratio = rupt_aspect_ratio
		self.tectonic_region_type = tectonic_region_type
		self.time_span = time_span
		self.grid_spacing = grid_spacing
		self.mesh_spacing = mesh_spacing
		self.area_grid = self.__create_grid()
		self.normalization_factor = self.__get_normalization_factor()
		
	def get_num_ruptures(self):
		num_grid_points = len(self.area_grid)
		num_mags = len(self.freq_mag_dist.getAnnualOccurrenceRates())
		num_nodal_planes = len(self.nodal_plane_pmf)
		num_hypo_depths = len(self.hypo_depth_pmf)
		return num_grid_points * num_mags * num_nodal_planes * num_hypo_depths
		
	def get_rupture(self,rupt_index):
		"""
		Create point source and return corresponding rupture.
		"""
		point_index, point_rupt_index = self.__get_point_rupture_indexes(rupt_index)
		
		point_weight = cos(radians(self.area_grid[point_index][1])) / self.normalization_factor
		
		point_source = PoissonianPointSource(Point(self.area_grid[point_index][0],self.area_grid[point_index][1]),
													self.freq_mag_dist.getWeightedMfd(point_weight),
													self.nodal_plane_pmf,
													self.hypo_depth_pmf,
													self.upper_seismo_depth,
													self.lower_seismo_depth,
													self.mag_scaling_rel,
													self.rupt_aspect_ratio,
													self.tectonic_region_type,
													self.time_span,
													self.mesh_spacing)
													
		return point_source.getRupture(point_rupt_index)
		
	def __get_normalization_factor(self):
		f = 0.0
		for lon, lat in self.area_grid:
			f = f + cos(radians(lat))
		return f
		
	def __get_point_rupture_indexes(self,rupt_index):
		"""
		Returns grid point index and corresponding rupture index.
		"""
		num_grid_points = len(self.area_grid)
		num_mags = len(self.freq_mag_dist.getAnnualOccurrenceRates())
		num_nodal_planes = len(self.nodal_plane_pmf)
		num_hypo_depths = len(self.hypo_depth_pmf)

		point_index = rupt_index % num_grid_points
		point_rupt_index = (rupt_index / num_grid_points) % (num_mags * num_nodal_planes * num_hypo_depths)

		return point_index, point_rupt_index
		
	def __create_grid(self):
		"""
		Create uniform grid over region.
		"""
		polygon_coords = self.__resample_polygon()
		
		# define polygon from polygon coords
		polygon = Polygon(polygon_coords)
		
		# get bounding box coordinates
		bounding_box = polygon.bounds
		min_lon = bounding_box[0]
		max_lon = bounding_box[2]
		min_lat = bounding_box[1]
		max_lat = bounding_box[3]
		
		# compute number of nodes along lat and lon
		n_lat = int(round((max_lat - min_lat) / self.grid_spacing)) + 1
		n_lon = int(round((max_lon - min_lon) / self.grid_spacing)) + 1
		
		# create grid of points inside polygon
		grid = []
		for i in range(n_lat):
			for j in range(n_lon):
				lat = min_lat + i * self.grid_spacing
				lon = min_lon + j * self.grid_spacing
				p = shapelyPoint(lon,lat)
				if polygon.contains(p) or polygon.touches(p):
					grid.append((lon,lat))
		return grid
		
	def __resample_polygon(self):
		"""
		Resample polygon with points that are equally spaced in each polygon segment.
		The spacing is approximately equal to GREAT_CIRCLE_DELTA. The spacing is not
		exactly equal to GREAT_CIRCLE_DELTA because each polygon segment is not exactly a
		multiple of the GREAT_CIRCLE_DELTA. 
		NOTE: The method assumes the first and last points in the polygon boundary
		to be same. (check for that)
		"""
		# loop over polygons segments. Resample each segment
		# in n equally spaced points along the geodesics connecting
		# the first and last points of the polygon segment.
		# n is computed as:
		# n = int(round(segment_length / GREAT_CIRCLE_DELTA))+1
		resampled_polygon = [self.area_boundary[0]]
		for i in range(len(self.area_boundary) - 1):
			p1 = Point(self.area_boundary[i][0],self.area_boundary[i][1])
			p2 = Point(self.area_boundary[i+1][0],self.area_boundary[i+1][1])
			segment_length = p1.getHorizontalDistance(p2)
			n = int(round(segment_length / GREAT_CIRCLE_DELTA))+1
			if n > 3:
				# -2 because the method computes only the intemediate points,
				# excluding the first and last points
				intermediate_points = g.npts(p1.longitude,p1.latitude,p2.longitude,p2.latitude,n - 2)
				resampled_polygon.extend(intermediate_points)
			resampled_polygon.append(self.area_boundary[i+1])
		return resampled_polygon