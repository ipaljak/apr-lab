from chromosome import *
from random import shuffle

class Genetic:
	def __init__(self, f, n, l, h, p_m, cross, p, gen, t, pr):
		self.f = f
		self.n = n
		self.l = l
		self.h = h
		self.p_m = p_m
		self.cross = cross
		self.p = p
		self.gen = gen
		self.population = []
		self.t = t
		self.pr = pr

	def initialize(self):
		raise NotImplementedError()		

	def work(self):
		self.initialize()

		for i in range(self.p):
			self.population[i].calculate_fitness(self.f)

		for i in range(self.gen):
			best = 0
			for j in range(self.p):
				if self.population[j].fitness < self.population[best].fitness:
					best = j

			if i == 0 or self.population[best].fitness < past_best:
				if self.pr:
					print "Generation:", i, "| Value:", self.population[best].fitness
					print self.population[best].values
				past_best = self.population[best].fitness
			
			self.get_next_generation()

		return self.population[best]

	def get_next_generation(self):
		ind = list(range(self.p))
		shuffle(ind)
		ind = ind[0:self.t]
		ind.sort(key=lambda x : self.population[x].fitness)

		child = self.cross.get_child(self.population[ind[0]], self.population[ind[1]])
		child.mutate(self.p_m)
		child.calculate_fitness(self.f)
		self.population[ind[self.t - 1]] = child

class GeneticDouble(Genetic):
	def __init__(self, f, n, l, h, p_m, cross, p, gen, t, pr=True):
		Genetic.__init__(self, f, n, l, h, p_m, cross, p, gen, t, pr)

	def initialize(self):
		self.population = []
		for i in range(self.p):
			self.population.append(ChromosomeDouble(self.n, self.l, self.h))

class GeneticBinary(Genetic):
	def __init__(self, f, n, l, h, p_m, cross, p, gen, t, m, pr=True):
		Genetic.__init__(self, f, n, l, h, p_m, cross, p, gen, t, pr)
		self.m = m

	def initialize(self):
		self.population = []
		for i in range(self.p):
			self.population.append(ChromosomeBinary(self.n, self.l, self.h, self.m))			

class CrossDouble1:
	def get_child(self, p1, p2):
		return ChromosomeDouble(p1.n, p1.l, p1.h, [(v1 + v2) / 2 for v1, v2 in zip(p1.values, p2.values)])

class CrossDouble2:
	def get_child(self, p1, p2):
		v = []
		for i in range(p1.n):
			if i % 2 == 0:
				v.append(p1.values[i])
			else:
				v.append(p2.values[i])
		return ChromosomeDouble(p1.n, p1.l, p1.h, v)

class CrossBinary1:
	def get_child(self, p1, p2):
		return ChromosomeBinary(p1.n, p1.l, p1.h, p1.m, [v1 & v2 for v1, v2 in zip(p1.values, p2.values)])		

class CrossBinary2:
	def get_child(self, p1, p2):
		return ChromosomeBinary(p1.n, p1.l, p1.h, p1.m, [v1 | v2 for v1, v2 in zip(p1.values, p2.values)])
		
