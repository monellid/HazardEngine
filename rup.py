from pyproj import Geod
import numpy

class ProbEqkRupture:
	"""
	Class defining probabilistic earthquake rupture.
	"""
	
	def __init__(self,magnitude,rupture_surface,hypocenter,rake,tectonic_region_type,probability_occurrence):
		"""
		magnitude: earthquake magnitude
		rupture_surface: 2D array representing rupture surface mesh
		hypocenter: hypocenter location (longitude,latitude,depth)
		rake: rake angle
		tectonic_region_type: tectonic region type
		probability_occurrence: probability of occurrence 
		"""
		self.magnitude = magnitude
		self.rupture_surface = rupture_surface
		self.hypocenter = hypocenter
		self.rake = rake
		self.tectonic_region_type = tectonic_region_type
		self.probability_occurrence = probability_occurrence
		
	def getShortestDistance(self,point):
		"""
		Compute shortest distance (in km) between point and rupture.
		The shortest distance is defined as the minimum
		distance between the given point and the points
		constituting the rupture surface mesh.
		"""
		g = Geod(ellps="sphere")
		point_list = self.rupture_surface.ravel()
		lons_rup = numpy.array([point_list[i].longitude for i in range(len(point_list))])
		lats_rup = numpy.array([point_list[i].latitude for i in range(len(point_list))])
		lons_point = numpy.array([point.longitude for i in range(len(point_list))])
		lats_point = numpy.array([point.latitude for i in range(len(point_list))])
		fwd_azs,back_azs,hor_dists = g.inv(lons_point,lats_point,lons_rup,lats_rup)
		vert_dists = numpy.array([(point.depth - point_list[i].depth) * 1e3 for i in range(len(point_list))])
		dists = numpy.sqrt(hor_dists**2 + vert_dists**2)
		return min(dists) * 1e-3
		
	def getJoynerBooreDistance(self,point):
		"""
		Compute shortest distance (in km) between point and
		rupture surface projection. The shortest distance is
		defined as the minimum distance between the point and
		and the surface projections of the points constituing
		the rupture surface mesh.
		Note that if the point lies inside the surface 
		projection of the rupture boundary, the distance in zero.
		"""
		
	def getDistanceX(self,point):
		"""
		The shortest horizontal distance from a Site to a line 
		defined by extending the top edge of the rupture to 
		infinity in both directions. Values on the 
		hanging-wall are positive and those on the foot-wall 
		are negative.
		"""