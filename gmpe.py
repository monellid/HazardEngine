from math import *

class SadighEtAl1997:
	"""
	Class implementing Sadigh Et. Al. 1997 GMPE.
	For the moment only PGA, on rock, for the
	average horizontal component. This is what
	is needed by the PEER Tests.
	"""
	
	def __init__(self,std_type,trucation_level,max_distance):
		"""
		std_type: standard deviation type (None, Total)
		truncation_level: number of standard deviations for truncation
		max_distance: maximum distance beyond which median
		ground motion value is zero (mean -> -infinity)
		"""
		self.std_type = std_type
		self.truncation_level = truncation_level
		self.max_distance = max_distance
	
	def getMean(self,rupture,point):
		"""
		Returns mean PGA.
		"""
		mag = rupture.magnitude
		dist = rupture.getShortestDistance(point)
		if rupture.rake > 45 and rupture.rake < 135:
			faulting_style = 'REVERSE'
		else:
			faulting_style = 'OTHER'
			
		if dist > max_distance:
			return -float('Inf')
		else:
			c1_rlt
			c2_rlt
			c3
			c4
			c5_rlt
			c6_rlt
			c7_r
			if mag <= 6.5:
				mean = coeff.c1_rlt + c2_rlt * mag +
				coeff.c3 * (Math.pow( (8.5 - mag), 2.5)) +
				coeff.c4 * (Math.log(dist + Math.exp(c5_rlt + c6_rlt * mag))) +
				coeff.c7_r * (Math.log(dist + 2));
			else:
				mean = coeff.c1_rgt + c2_rgt * mag +
				coeff.c3 * (Math.pow( (8.5 - mag), 2.5)) +
				coeff.c4 * (Math.log(dist + Math.exp(c5_rgt + c6_rgt * mag))) +
				coeff.c7_r * (Math.log(dist + 2));