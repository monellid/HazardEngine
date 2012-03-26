from pyproj import Geod
from math import *
import numpy


class Point:
	"""
	Class representing a geographical point in terms of longitude, latitude, both in decimal degrees,
	and depth (with respect to Earth surface, in km).
	Depth = 0 identifies a point on the Earth surface.
	Depth > 0 identifies a point below the Earth surface.
	Depth < 0 identifies a point above the Earth surface.
	"""
	
	#: The distance between two points for them to be considered equal,
    #: in km.
	EQUALITY_DISTANCE = 1e-3

	def __init__(self,longitude,latitude,depth=0):
		"""
		Create point from longitude, latitude, both in decimal degrees, and depth values (in km).
		"""
		self.longitude = longitude
		self.latitude = latitude
		self.depth = depth
		
	def __eq__(self, other):
		if other == None:
			return False
		return abs(self.getDistance(other)) <= self.EQUALITY_DISTANCE

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
		if total_distance == 0:
			return points
					
		horizontal_distance = self.getHorizontalDistance(point)
		azimuth = self.getAzimuth(point)
		
		bearing_angle = asin(horizontal_distance / total_distance)
		if point.depth != self.depth:
			sign = (point.depth - self.depth) / fabs(point.depth - self.depth) # if positive -> pointing downwards; if negative -> pointing upwards
		else:
			sign = 1
		vertical_increment = sign * distance * cos(bearing_angle)
		horizontal_increment = distance * sin(bearing_angle)

		number_locations = int(round(total_distance / distance) + 1)
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
		
	def getAverageAzimuth(self):
		"""
		The average azimuth is computed following the approach
		described in http://en.wikipedia.org/wiki/Mean_of_circular_quantities.
		That is segment's azimuths and lengths are interpreted as
		polar coordinates and converted to cartesian coordinates
		(therefore defining a set of cartesian vectors). The mean azimuth is given
		by the angle of the vector obtained by summing all vectors.
		"""
		#TODO: if line is perfectly vertical return zero
		
		if len(self.point_list) == 2:
			return self.point_list[0].getAzimuth(self.point_list[1])
		else:
			# loop over line segments and compute
			# segment's azimuths and lenghts
			azimuths = []
			lengths = []
			for i in range(len(self.point_list) - 1):
				azimuths.append(radians(self.point_list[i].getAzimuth(self.point_list[i+1])))
				lengths.append(self.point_list[i].getHorizontalDistance(self.point_list[i+1]))
			total_length = sum(lengths)				
			
			# convert from polar to cartesian coordinates
			vectors = []
			for i in range(len(azimuths)):
				vectors.append(self.__getCartesianVector(lengths[i],azimuths[i]))
				
			# sum all vectors. this represents the mean direction,
			# from which we can extract the mean angle
			v = vectors[0]
			for i in range(1,len(vectors)):
				v = v + vectors[i]
				
			# extract angle
			strike = degrees(atan(v[0] / v[1]))
			
			if strike < 0:
				strike = strike + 360.0
				
			if strike >= 360.0:
				strike = strike - 360.0
			
			return strike
			
	def __getCartesianVector(self,radius,azimuth):
		"""
		Return cartesian vector from polar vector.
		Azimuth is measured from the y axis, and
		is assumed to be in radians.
		"""
		x = radius * sin(azimuth)
		y = radius * cos(azimuth)

		return numpy.array([x,y])

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
		
	def getLength(self):
		"""
		Computes line length.
		"""
		length = 0.0
		for i in range(len(self.point_list) - 1):
			dist = self.point_list[i].getDistance(self.point_list[i+1])
			length = length + dist
		
		return length
		
	def getResampledLineInNpoints(self,n_points):
		"""
		Resample line into n points.
		Points in the resampled line are equally spaced
		along the original line.
		"""
		resampled_line = [self.point_list[0]]
		
		if n_points == 1:
			return Line(resampled_line)
		
		# compute section length compatible with n_points
		line_length = self.getLength()
		section_length = line_length / (n_points - 1)
		
		for i in range(n_points - 1):
			tot_length = (i+1) * section_length
			p1,p2,offset = self.getClosestPointsAndOffset(tot_length)
			azimuth = p1.getAzimuth(p2)
			horizontal_distance = p1.getHorizontalDistance(p2)
			vertical_distance = p2.depth - p1.depth
			sign = 1
			if p2.depth - p1.depth < 0.0:
				sign = -1
			total_distance = sqrt(horizontal_distance**2 + vertical_distance**2)
			# bearing angle: angle between line connecting p1 and p2
			# and great circle line passing through p1 and with same azimuth
			bearing_angle = acos(horizontal_distance / total_distance)
			horizontal_distance = offset * cos(bearing_angle)
			vertical_distance = offset * sin(bearing_angle) * sign
			p = p1.getPoint(horizontal_distance,vertical_distance,azimuth)
			resampled_line.append(p)
		
		return Line(resampled_line)
			
	def getClosestPointsAndOffset(self,length):
		"""
		return the two points in a line
		which define lengths (from the first point of the line)
		that are the closest to 'length'.
		Returns also the offset between lenght and length as
		given by the first point.
		"""
		#TODO: check length >= 0
		if length >= self.getLength():
			return self.point_list[-2],self.point_list[-1],self.point_list[-2].getDistance(self.point_list[-1])	
		elif length == 0:
			return self.point_list[0],self.point_list[1],0.0 
		else:
			for i in range(len(self.point_list) - 1):
				length1 = Line(self.point_list[:i+1]).getLength()
				length2 = Line(self.point_list[:i+2]).getLength()
				if length >= length1 and length < length2:
					return self.point_list[i],self.point_list[i+1],length - length1
					
	def getClosestPoint(self,length):
		"""
		Returns the point in a line
		that corresponds to a length (from the first point of the line)
		that is closest to 'length'
		"""
		p1,p2,offset = self.getClosestPointsAndOffset(length)
		dist = p1.getDistance(p2)
		
		if offset >= dist:
			return p2
		else:
			closest_p = p1
			if dist - offset < offset:
				closest_p = p2
			return closest_p
					
class Triangle:
	"""
	Class defining a triangle, as defined by three geographical locations.
	"""
	def __init__(self,p1,p2,p3):
		"""
		Create triangle from three points.
		"""
		# TODO: check the three points are all different
		self.p1 = p1
		self.p2 = p2
		self.p3 = p3
		
	def getArea(self):
		"""
		Returns triangle area (in km**2) using Heron's formula
		http://mathworld.wolfram.com/HeronsFormula.html
		"""
		# get position vectors for the three points
		P1 = getPositionVector(self.p1.longitude,self.p1.latitude,self.p1.depth)
		P2 = getPositionVector(self.p2.longitude,self.p2.latitude,self.p2.depth)
		P3 = getPositionVector(self.p3.longitude,self.p3.latitude,self.p3.depth)
		
		# get vectors joing triangle vertices
		P1P2 = P2 - P1
		P2P3 = P3 - P2
		P3P1 = P1 - P3
		
		# compute norms (that is sides lengths)
		p1p2 = sqrt(numpy.dot(P1P2,P1P2))
		p2p3 = sqrt(numpy.dot(P2P3,P2P3))
		p3p1 = sqrt(numpy.dot(P3P1,P3P1))
		
		# compute semiperimeter
		s = (p1p2 + p2p3 + p3p1) / 2
		
		# compute area
		a = sqrt(s * (s - p1p2) * (s - p2p3) * (s - p3p1))
		
		return a
		
	def getInclination(self):
		"""
		Returns angle between Earth surface and plane
		containing the triangle (between 0 and 90 degrees).
		"""
		# get position vectors for the three points
		P1 = getPositionVector(self.p1.longitude,self.p1.latitude,self.p1.depth)
		P2 = getPositionVector(self.p2.longitude,self.p2.latitude,self.p2.depth)
		P3 = getPositionVector(self.p3.longitude,self.p3.latitude,self.p3.depth)
		
		# get two vectors joing one vertex
		# (p1) with the remaining two (p2, p3)
		P1P2 = P2 - P1
		P1P3 = P3 - P1
		
		# make the cross product of the two vectors,
		# the result is a vector perpendicular to the
		# triangle
		normal = numpy.cross(P1P2,P1P3)
		normal = normal / sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
		# compute unit vector normal to earth surface at p1
		normal_p1 = P1 / sqrt(P1[0]**2 + P1[1]**2 + P1[2]**2)
		# compute angle in between: this is the dip
		dip = degrees(acos(numpy.dot(normal_p1,normal)))
		
		# this is to take into account the case in which
		# the vector perpendicular to the triangle is
		# not pointing towards the Earth surface, but in the
		# opposite direction.
		if dip > 90.0:
			dip = 180 - dip
			
		return dip
		
	def getCentroid(self):
		"""
		Returns triangle centroid (as cartesian vector).
		Implements equation (5) in 
		http://mathworld.wolfram.com/GeometricCentroid.html 
		"""
		# get position vectors for the three points
		P1 = getPositionVector(self.p1.longitude,self.p1.latitude,self.p1.depth)
		P2 = getPositionVector(self.p2.longitude,self.p2.latitude,self.p2.depth)
		P3 = getPositionVector(self.p3.longitude,self.p3.latitude,self.p3.depth)
		
		return (P1+P2+P3) / 3.0
		
class Rectangle:
	"""
	Class defining a rectangular plane
	in 3D as defined by four geographical locations.
	"""
	def __init__(self,p1,p2,p3,p4):
		"""
		p1 = upper left corner
		p2 = upper right corner
		p3 = lower left corner
		p4 = lower right corner
		"""
		self.p1 = p1
		self.p2 = p2
		self.p3 = p3
		self.p4 = p4
		
	def getMinimumDistance(self,p):
		"""
		Computes minimum distance (in km) between point (p)
		and rectangle.
		The algorithm works as follows:
		1) compute point projection on the plane
		containing the rectangle
		2) check if the point lies in or outside
		the rectangle
		3) compute minimum distance:
			if inside the rectangle compute minimum distance as:
				equation 13 in http://mathworld.wolfram.com/Point-PlaneDistance.html
			if outside compute minimum distance in the following way:
				i) find closest segment
				ii) compute distance between point and closest segment
		"""
		
		# find plane equation coefficients (equations 19-20 in http://mathworld.wolfram.com/Plane.html)
		P1 = getPositionVector(self.p1.longitude,self.p1.latitude,self.p1.depth)
		P2 = getPositionVector(self.p2.longitude,self.p2.latitude,self.p2.depth)
		P3 = getPositionVector(self.p3.longitude,self.p3.latitude,self.p3.depth)
		P4 = getPositionVector(self.p4.longitude,self.p4.latitude,self.p4.depth)
		M = numpy.array([P1,P2,P3])
		M1 = numpy.array(M)
		M2 = numpy.array(M)
		M3 = numpy.array(M)
		M1[:,0] = numpy.ones(3)
		M2[:,1] = numpy.ones(3)
		M3[:,2] = numpy.ones(3)
		A1 = numpy.linalg.det(M1)
		A2 = numpy.linalg.det(M2)
		A3 = numpy.linalg.det(M3)
		A = numpy.linalg.det(M)
		C = numpy.array([A1,A2,A3])
		
		# find projection of p onto plane containing the rectangle
		P = getPositionVector(p.longitude,p.latitude,p.depth)
		factor = (A - numpy.sum(C * P)) / numpy.sum(C**2)
		P0 = P + C * factor
		
		# translate coordinate system origin (at the Earth center)
		# to the rupture bottom left corner (that is P3)
		# n 
		P1n = P1 - P3
		P2n = P2 - P3
		P4n = P4 - P3
		P0n = P0 - P3
		P3n = P3 - P3
		
		# rotate coordinate system to have x axis aligned along rupture bottom edge
		# and y axis aligned along rupture left edge
		# for the equation see http://mathworld.wolfram.com/EulerAngles.html
		vec1 = P4n / numpy.sqrt(numpy.sum(P4n * P4n))
		vec1[2] = 0.0
		vec2 = P1n / numpy.sqrt(numpy.sum(P1n * P1n))
		vec2[2] = 0.0
		phi = numpy.arccos(numpy.dot(numpy.array([1.0,0.0,0.0]),vec1))
		theta = numpy.arccos(numpy.dot(P1n/numpy.sqrt(numpy.sum(P1n * P1n)),vec2))
		psi = numpy.arccos(numpy.dot(vec1,P4n/numpy.sqrt(numpy.sum(P4n * P4n))))
		rot = numpy.array([[cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi),cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi),sin(psi)*sin(theta)],
							[-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi),-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi),cos(psi)*sin(theta)],
							[sin(theta)*sin(phi),-sin(theta)*cos(phi),cos(theta)]])
		P1n = numpy.dot(rot,P1n)
		P2n = numpy.dot(rot,P2n)
		P4n = numpy.dot(rot,P4n)
		P0n = numpy.dot(rot,P0n)
		
		# check if point is inside or outside rectangle
		inside = False
		if P0n[0] >= P3n[0] and P0n[0] <= P4n[0] and P0n[1] >= P3n[1] and P0n[1] <= P1n[1]:
			inside = True
			
		# compute minimum distance
		if inside:
			return numpy.sqrt(numpy.sum((P0 - P)**2))
		else:
			
			# compute distance based on closest segment/point
			if P0n[0] <= P1n[0] and P0n[1] >= P1n[1]:
				# closest point is P1
				return numpy.sqrt(numpy.sum((P0n - P1n)**2)) + numpy.sqrt(numpy.sum((P0 - P)**2))
			
		
def getPositionVector(longitude,latitude,depth):
	"""
	Returns the position vector (in Cartesian coordinates,km) of
	a geographical location defined by longitude, latitude and
	depth.
	For the equation see:
	http://mathworld.wolfram.com/SphericalCoordinates.html
	Longitude and latitudes are supposed to be in degrees and 
	depth in km."""
	# TODO: check that depth is < than earth radius. Earth
	# radius should be defined in a common place.
	earth_radius = 6371 # in km
	theta = radians(longitude)
	phi = radians(90.0 - latitude)
	x = (earth_radius - depth) * cos(theta) * sin(phi)
	y = (earth_radius - depth) * sin(theta) * sin(phi)
	z = (earth_radius - depth) * cos(phi)
	return numpy.array([x,y,z])
	
def getSphericalPositionVector(x,y,z):
	"""
	Returns the position vector (in Spherical coordinates,
	latitude, longitude in degrees, and depth in km).
	For the equation see:
	http://mathworld.wolfram.com/SphericalCoordinates.html
	x,y,z are supposed to be in km.
	"""
	earth_radius = 6371 # in km
	numpy.seterr(divide='ignore')
	radius = sqrt(x**2 + y**2 + z**2) # TODO: check that radius is not zero, in this case the spherical vector is not defined
	theta = atan(numpy.divide(y,x))
	phi = acos(z / radius) 
	
	longitude = degrees(theta)
	latitude = 90.0 - degrees(phi)
	depth = earth_radius - radius
	
	return numpy.array([longitude,latitude,depth])