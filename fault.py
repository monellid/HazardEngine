from math import *
from geo import *
import numpy
		
class SimpleFaultSurface:
	"""
	Class defining Simple Fault Surface.
	"""
	
	def __init__(self,fault_trace,upper_seismo_depth,lower_seismo_depth,dip,mesh_spacing = 1):
		"""
		Represents fault surface as regular 3D mesh of points derived from:
		fault_trace: list of points representing intersection between fault surface and Earth surface
		upper_seismo_depth: upper seismogenic depth (i.e. depth to fault's top edge, in km)
		lower_seismo_depth: lower seismogenic depth (i.e. depth to fault's bottom edge, in km)
		dip: dip angle (in degrees)
		mesh_spacing: spacing (in km) between mesh points
		"""
		# TODO: check that all points in fault_trace have depth = 0
		# TODO: check upper_seismo_depth >= lower_seismo_depth >= 0
		# TODO: check dip >0 && dip <=90
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
		top_edge = self.getFaultTopEdge()
		
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
		
	def getFaultTopEdge(self):
		"""
		Computes fault top edge coordinates by translating the fault trace along dip
		from the Earth surface to the upper seismogenic depth in a direction 
		perpendicular to the fault strike (computed considering the first and last points).
		Point's coordinates are then resampled to be equally spaced (with spacing equal
		to mesh_spacing).
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
			
		top_edge = Line(top_edge).getResampledLine(self.mesh_spacing).point_list
		top_edge.append(top_edge[-1].getPoint(self.mesh_spacing,0.0,top_edge[-2].getAzimuth(top_edge[-1]))) # this to get exact match with opensha
		return top_edge#Line(top_edge).getResampledLine(self.mesh_spacing).point_list
		
class PoissonianFaultSource:
	"""
	Class defining Poissonian Fault source, that is a fault source with
	ruptures having Poissonian probability of occurrence.
	"""
	
	def __init__(self,fault_surf,freq_mag_dist,mag_scaling_rel,rake,rup_aspect_ratio):
		"""
		fault_surf: fault surface
		freq_mag_dist: frequency magnitude distribution
		mag_scaling_rel: magnitude scaling relationship
		rake: rake angle
		rup_aspect_ratio: rupture aspect ratio
		"""
		self.fault_surf = fault_surf
		self.freq_mag_dist = freq_mag_dist
		self.mag_scaling_rel = mag_scaling_rel
		self.rake = rake
		self.rup_aspect_ratio = rup_aspect_ratio
		
	def getNumRuptures(self):
		"""
		Count number of ruptures. Loop over magnitude values in the
		frequency magnitude distribution. For each magnitude, computes
		rupture's dimensions (length and width). Based on rupture's
		dimensions and rupture offset (assumed equal to surface mesh spacing), 
		computes how many ruptures can be placed on the fault surface. 
		"""
		num_ruptures = 0.0
		
		fault_length = (self.fault_surf.surface.shape[1] - 1) * self.fault_surf.mesh_spacing
		fault_width = (self.fault_surf.surface.shape[0] - 1) * self.fault_surf.mesh_spacing
		num_mesh_points = self.fault_surf.surface.shape[0] * self.fault_surf.surface.shape[1]
		
		occurrence_rates = self.freq_mag_dist.getAnnualOccurrenceRates()
		for mag,rate in occurrence_rates:
			
			# compute rupture dimensions
			area = self.mag_scaling_rel.getMedianArea(mag)
			rup_length = sqrt(area * self.rup_aspect_ratio)
			rup_width = area / rup_length
			
			# round rupture dimensions with respect to mesh_spacing
			rup_length = round(rup_length / self.fault_surf.mesh_spacing) * self.fault_surf.mesh_spacing
			rup_width = round(rup_width / self.fault_surf.mesh_spacing) * self.fault_surf.mesh_spacing
			
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
			# dimensions are greater than fault dimensions
			if rup_length > fault_length and rup_width > fault_width:
				rup_length = fault_length
				rup_width = fault_width
			
			if rup_length < self.fault_surf.mesh_spacing and rup_width < self.fault_surf.mesh_spacing:
				num_ruptures = num_ruptures + num_mesh_points
			elif rup_length >= fault_length and rup_width >= fault_width:
				num_ruptures = num_ruptures + 1
			else:
				# computes number of ruptures along length and along width
				num_rup_along_length = int((fault_length - rup_length) / self.fault_surf.mesh_spacing + 1)
				num_rup_along_width = int((fault_width - rup_width) / self.fault_surf.mesh_spacing + 1)
				num_ruptures = num_ruptures + num_rup_along_length * num_rup_along_width
		
		return num_ruptures