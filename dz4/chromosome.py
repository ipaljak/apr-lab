from random import random, randint

class Chromosome:
	def __init__(self, n, l, h):
		self.n = n
		self.l = l
		self.h = h
		self.fitness = 0

	def calculate_fitness(self, f):
		self.fitness = self.evaluate(f)

	def __lt__(self, other):
		return self.fitness < other.fitness

class ChromosomeDouble(Chromosome):
		def __init__(self, n, l, h, values=[]):
			Chromosome.__init__(self, n, l, h)
			if len(values) > 0:
				self.values = values
				return
			self.values = []
			for i in range(n):
				self.values.append(self.l + random() * (self.h - self.l))

		def mutate(self, p):
			for i in range(self.n):
				if random() <= p:
					self.values[i] = self.l + random() * (self.h - self.l)

		def evaluate(self, f):
			return f.evaluate(self.values)					

class ChromosomeBinary(Chromosome):
		def __init__(self, n, l, h, m, values=[]):
			Chromosome.__init__(self, n, l, h)
			self.m = m
			if len(values) > 0:
				self.values = values
				return
			self.values = []
			for i in range(n):
				self.values.append(randint(0, 2**self.m - 1))

		def mutate(self, p):
			for i in range(self.n):
				for j in range(self.m):
					if random() <= p:
						self.values[i] ^= 1 << j

		def evaluate(self, f):
			return f.evaluate([self.l + float(x) / (2**self.m - 1) * (self.h - self.l) for x in self.values])
