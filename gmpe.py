from math import *

class SadighEtAl1997:
	"""
	Class implementing Sadigh Et. Al. 1997 GMPE.
	"""
	
	def __init__(self,std_type,trucation_level):
		"""
		std_type: standard deviation type (None, Total)
		truncation_level: number of standard deviations for truncation
		"""
		self.std_type = std_type
		self.truncation_level = truncation_level
	
	def getMean(self,rupture,site):
		"""
		Returns mean PGA.
		"""
		
		