from pyproj import Geod
from math import *

class Point:
	"""
	Class representing a geographical point in terms of longitude, latitude, both in decimal degrees,
	and depth (with respect to Earth surface, in km).
	Depth = 0 identifies a point on the Earth surface.
	Depth > 0 identifies a point below the Earth surface.
	Depth < 0 identifies a point above the Earth surface.
	"""

	def __init__(self,longitude,latitude,depth=0):
		"""
		Create point from longitude, latitude, both in decimal degrees, and depth values (in km).
		"""
		self.longitude = longitude
		self.latitude = latitude
		self.depth = depth

	def getPoint(self,horizontal_distance,vertical_distance,azimuth):
		"""
		Get point with given horizontal, and vertical distances (in km, 
		vertical distance: positive-downward, negative-upward) 
		and azimuth (in degrees) from current point. 
		"""
		# TODO: check horizontal distance is positive
		g = Geod(ellps="sphere")
		longitude, latitude, back_azimuth = g.fwd(self.longitude, self.latitude, azimuth, horizontal_distance * 1e3) # 1e3 is needed to convert from km to m
		depth = self.depth + vertical_distance
		return Point(longitude,latitude,depth)
		
	def getAzimuth(self,point):
		"""
		Get azimuth (in degrees) between current point and provided point (point).
		"""
		g = Geod(ellps="sphere")
		forward_azimuth,back_azimuth,distance = g.inv(self.longitude, self.latitude,point.longitude, point.latitude)
		return forward_azimuth
		
	def getHorizontalDistance(self,point):
		"""
		Get horizontal distance (great circle distance, in km) between current point
		and provided point (point).
		"""
		g = Geod(ellps="sphere")
		forward_azimuth,back_azimuth,horizontal_distance = g.inv(self.longitude, self.latitude,point.longitude, point.latitude)
		return horizontal_distance * 1e-3 # 1e-3 is needed to convert from m to km
		
	def getDistance(self,point):
		"""
		Get distance (in km) between current point and provided point (point).
		Distance is calculated using pythagoras theorem, where the hypotenuse is
		the distance and the other two sides are the horizontal distance
		(great circle distance) and vertical distance (depth difference between 
		the two locations).
		"""
		horizontal_distance = self.getHorizontalDistance(point)
		vertical_distance = point.depth - self.depth
		return sqrt(horizontal_distance**2 + vertical_distance**2)
		
	def getEquallySpacedPoints(self,point,distance):
		"""
		Get list of points equally spaced (spacing equal to distance [in km])
		between current point and provided point (point).
		"""
		points = []
		points.append(self)
		
		total_distance = self.getDistance(point)		
		horizontal_distance = self.getHorizontalDistance(point)
		azimuth = self.getAzimuth(point)
		
		bearing_angle = asin(horizontal_distance / total_distance)
		if point.depth != self.depth:
			sign = (point.depth - self.depth) / fabs(point.depth - self.depth) # if positive -> pointing downwards; if negative -> pointing upwards
		else:
			sign = 1
		vertical_increment = sign * distance * cos(bearing_angle)
		horizontal_increment = distance * sin(bearing_angle)
		number_locations = int(total_distance / distance) + 1
		for i in range(1,number_locations):
			p = points[-1]
			points.append(p.getPoint(horizontal_increment,vertical_increment,azimuth))
		return points
		
class Line:
	"""
	Class defining a geographical line (that is a list of geographical points)
	"""

	def __init__(self,point_list):
		"""
		Defines line as list of points.
		"""
		# TODO: check that there are at least two points and they are not coincident
		# TODO: check that the line does not intersect itself
		self.point_list = point_list

	def getResampledLine(self,section_length):
		"""
		Resamples line into sections with length equal to 
		section_length.
		"""
		resampled_line = []
		# First resample first section. Then loop over remaining points in the fault trace
		# and resample remaining sections. Extend the list with the resampled points,
		# except the first (because it's already contained in the previous set of 
		# resampled points)
		resampled_line.extend(self.point_list[0].getEquallySpacedPoints(self.point_list[1],section_length))
		num_points = len(self.point_list)
		for i in range(2,num_points):
			points = resampled_line[-1].getEquallySpacedPoints(self.point_list[i],section_length)
			resampled_line.extend(points[1:])
		return Line(resampled_line)