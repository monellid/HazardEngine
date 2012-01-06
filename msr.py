from math import *

class PeerTestMagAreaScalingRel:
	"""
	Class implementing Peer Tests magnitude-area scaling relationship.
	"""
	def getMedianArea(self,mag):
		"""
		Returns median area (in km**2) from magnitude (mag).
		"""
		return pow(10,mag - 4.0)