from math import sqrt, sin

class Function:
	def evaluate(self, x):
		raise NotImplementedError()

class F1(Function):
	def evaluate(self, x):
		return 100 * (x[1] - x[0]**2)**2 + (1 - x[0])**2

class F3(Function):
	def evaluate(self, x):
		return sum([(a - i - 1)**2 for i, a in enumerate(x)])

class F6(Function):
	def evaluate(self, x):
		s = sum([a**2 for a in x])
		return 0.5 + (sin(sqrt(s))**2 - 0.5) / (1 + 0.001 * s)**2

class F7(Function):
	def evaluate(self, x):
		s = sum([a**2 for a in x])
		return s**0.25 * (1 + sin(50 * s**0.1)**2)
