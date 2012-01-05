from pyproj import Geod
from math import *

class Location:
	"""
	Class representing geographical location in terms of longitude, latitude, both in decimal degrees,
	and depth (with respect to Earth surface, in km).
	Depth = 0 identifies a point on the Earth surface.
	Depth > 0 identifies a point below the Earth surface.
	Depth < 0 identifies a point above the Earth surface.
	"""

	def __init__(self,longitude,latitude,depth=0):
		"""
		Create Location from longitude, latitude, both in decimal degrees, and depth values (in km).
		"""
		self.longitude = longitude
		self.latitude = latitude
		self.depth = depth

	def getLocation(self,horizontal_distance,vertical_distance,azimuth):
		"""
		Get Location with given horizontal, and vertical distances (in km, 
		vertical distance: positive-downward, negative-upward) 
		and azimuth (in degrees) from current location. 
		"""
		# TODO: check horizontal distance is postive
		g = Geod(ellps="sphere")
		longitude, latitude, back_azimuth = g.fwd(self.longitude, self.latitude, azimuth, horizontal_distance * 1e3) # 1e3 is needed to convert from km to m
		depth = self.depth + vertical_distance
		return Location(longitude,latitude,depth)
		
	def getAzimuth(self,loc):
		"""
		Get azimuth (in degrees) between current location and provided location (loc).
		"""
		g = Geod(ellps="sphere")
		forward_azimuth,back_azimuth,distance = g.inv(self.longitude, self.latitude,loc.longitude, loc.latitude)
		return forward_azimuth
		
	def getHorizontalDistance(self,loc):
		"""
		Get horizontal distance (great circle distance, in km) between current location
		and provided location (loc).
		"""
		g = Geod(ellps="sphere")
		forward_azimuth,back_azimuth,horizontal_distance = g.inv(self.longitude, self.latitude,loc.longitude, loc.latitude)
		return horizontal_distance * 1e-3 # 1e-3 is needed to convert from m to km
		
		
	def getDistance(self,loc):
		"""
		Get distance (in km) between current location and provided location (loc).
		Distance is calculated using pythagoras theorem, where the hypotenuse is
		the distance and the other two sides are the horizontal distance
		(great circle distance) and vertical distance (depth difference between 
		the two locations).
		"""
		horizontal_distance = self.getHorizontalDistance(loc)
		vertical_distance = loc.depth - self.depth
		return sqrt(horizontal_distance**2 + vertical_distance**2)
		
	def getLocations(self,loc,distance):
		"""
		Get list of locations equally spaced (spacing equal to distance [in km])
		between current location and provided location (loc).
		"""
		locations = []
		locations.append(self)
		
		total_distance = self.getDistance(loc)		
		horizontal_distance = self.getHorizontalDistance(loc)
		angle = asin(horizontal_distance / total_distance)
		azimuth = self.getAzimuth(loc)
		if loc.depth != 0 or self.depth != 0:
			sign = (loc.depth - self.depth) / fabs(loc.depth - self.depth) # if positive -> pointing downwards; if negative -> pointing upwards
		else:
			sign = 1
		
		
		number_locations = int(total_distance / distance) + 1
		for i in range(1,number_locations):
			vertical_increment = sign * distance * cos(angle)
			horizontal_increment = distance * sin(angle)
			loc = locations[-1]
			locations.append(loc.getLocation(horizontal_increment,vertical_increment,azimuth))
		return locations
