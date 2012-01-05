from math import *
from location import Location

class FaultTrace:
	"""
	Class defining Fault Trace (2D line representing fault surface intersection with Earth surface).
	"""
	
	def __init__(self,location_list):
		"""
		Defines fault trace as list of locations (at the Earth surface).
		"""
		# TODO: check that there are at least two locations
		# TODO: check that all locations have depth = 0
		# TODO: check that that the fault trace line does not intersect itself
		self.location_list = location_list
		
	def getStrikeFromFirstAndLastPoints(self):
		"""
		Computes strike angle (in degrees) as the azimuth from the
		first and last points of the fault trace.
		"""
		return self.location_list[0].getAzimuth(self.location_list[-1])
		
	def getResampledTrace(self,section_length):
		"""
		Resample fault trace into sections with length equal to 
		section_length.
		"""
		num_points = len(self.location_list)
		for i in range(num_points-1):
			distance = self.location_list[i].getDistance(self.location_list[i+1])
			# compute number of sections
			num_sections
		
class SimpleFaultSurface:
	"""
	Class defining Simple Fault Surface (as a 3D mesh of locations).
	"""
	
	def __init__(self,fault_trace,upper_seismo_depth,lower_seismo_depth,dip):
		"""
		fault_trace: fault trace (FaultTrace)
		upper_seismo_depth: upper seismogenic depth (in km)
		lower_seismo_depth: lower seismogenic depth (in km)
		dip: dip angle (in degrees)
		"""
		# TODO: check upper_seismo_depth >= lower_seismo_depth >= 0
		# TODO: check dip >=0 && dip <=90
		self.fault_trace = fault_trace
		self.upper_seismo_depth = upper_seismo_depth
		self.lower_seismo_depth = lower_seismo_depth
		self.dip = dip
		
	def getFaultSurface(mesh_spacing):
		"""
		Compute 3D mesh representing fault surface. The mesh is
		constructed by translating the fault trace
		along the dip angle, on a direction perpendicular to
		the fault trace strike (computed considering the first
		and last locations).
		Distance between mesh points is given by the mesh_spacing
		parameter (in km).
		"""		
		# compute fault strike considering first and last points
		fault_strike = self.fault_trace.getStrikeFromFirstAndLastPoints()
		
		# compute fault top edge coordinates by translating the
		# fault trace from the Earth surface to the upper seismogenic depth
		# in a direction perpendicular to the fault strike
		fault_top_edge = []
		if self.dip < 90.0:
			horizontal_distance = self.upper_seismo_depth / atan(radians(self.dip))
		else:
			horizontal_distance = 0.0
		azimuth = fault_strike + 90.0
		for loc in self.fault_trace.location_list:
			fault_top_edge.append(loc.getLocation(horizontal_distance,self.upper_seismo_depth,azimuth))
			
		# resample fault trace into sections with length equal to mesh_spacing