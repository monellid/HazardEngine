from math import *
from geo import *
import numpy
		
class SimpleFaultSurface:
	"""
	Class defining Simple Fault Surface.
	"""
	
	def __init__(self,fault_trace,upper_seismo_depth,lower_seismo_depth,dip):
		"""
		fault_trace: list of points representing intersection between fault surface and Earth surface
		upper_seismo_depth: upper seismogenic depth (i.e. depth of fault's top edge, in km)
		lower_seismo_depth: lower seismogenic depth (i.e. depth of fault's bottom edge, in km)
		dip: dip angle (in degrees)
		"""
		# TODO: check that all points in fault_trace have depth = 0
		# TODO: check upper_seismo_depth >= lower_seismo_depth >= 0
		# TODO: check dip >0 && dip <=90
		self.fault_trace = fault_trace
		self.upper_seismo_depth = upper_seismo_depth
		self.lower_seismo_depth = lower_seismo_depth
		self.dip = dip
		
	def getFaultSurface(self,mesh_spacing = 1):
		"""
		Computes 3D mesh representing fault surface. The mesh is constructed by 
		translating the fault trace along the dip angle, along a direction perpendicular
		to the fault trace strike (computed considering the first and last points in the fault trace).
		Distance between mesh points is given by the mesh_spacing parameter (in km).
		Default value is 1 km (the error in reproducing the fault length and width
		is at most 1 km)
		"""		
		top_edge = self.getFaultTopEdge(mesh_spacing)
		
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
			points = point.getEquallySpacedPoints(bottom_point,mesh_spacing)
			mesh_points.extend(points)
		
		# organize mesh points into a 2D numpy array
		# number of rows corresponds to number of points along dip
		# number of columns corresponds to number of points along strike
		surface = numpy.array(mesh_points)
		surface = surface.reshape(len(top_edge),len(mesh_points)/len(top_edge))
		surface = numpy.transpose(surface)
		return surface
		
	def getFaultTopEdge(self,mesh_spacing):
		"""
		Computes fault top edge coordinates by translating the fault trace along dip
		from the Earth surface to the upper seismogenic depth in a direction 
		perpendicular to the fault strike (computed considering the first and last points).
		Points' coordinates are then recomputed so that points are equally spaced
		(with spacing equal to mesh_spacing).
		"""
		top_edge = []
		
		if self.dip < 90.0:
			horizontal_distance = self.upper_seismo_depth / atan(radians(self.dip))
		else:
			horizontal_distance = 0.0
		vertical_distance = self.upper_seismo_depth
		strike = self.fault_trace[0].getAzimuth(self.fault_trace[-1])
		azimuth = strike + 90.0
		for point in self.fault_trace:
			top_edge.append(point.getPoint(horizontal_distance,vertical_distance,azimuth))
		return Line(top_edge).getResampledLine(mesh_spacing).point_list