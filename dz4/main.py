from chromosome import *
from genetic import *
from function import *
import numpy
from sys import argv

def zad1():
	sol = []
	for f, n in zip([F1(), F3(), F6(), F7()], [2, 5, 2, 2]):
		g = GeneticDouble(f, n, -50, 150, 0.3, CrossDouble1(), 200, 100000, 3, 20)		
		s = g.work()
		sol.append(s.evaluate(f))
		print("----------------------------")
	for f, s in zip(["F1", "F3", "F6", "F7"], sol):
		print f, s

def zad2():
	sol = []
	for f in [F6(), F7()]:
		for n in [1, 3, 6, 10]:
			g = GeneticDouble(f, n, -50, 150, 0.2, CrossDouble1(), 200, 30000, 3)		
			s = g.work()
			sol.append((n, s.evaluate(f)))
	print "F6:"
	for n, s in sol[0:4]: print "   n =", n, "|", s
	print "F7:"
	for n, s in sol[4:8]: print "   n =", n, "|", s

def zad3():
	EPS = 10e-6
	for f in [F6(), F7()]:
		for n in [3, 6]:
			sol_d, sol_b = [], []
			for i in range(10):
				g = GeneticDouble(f, n, -50, 50, 0.2, CrossDouble1(), 200, 10000, 3, pr=False)
				s = g.work()
				sol_d.append(s.evaluate(f))
				g = GeneticBinary(f, n, -50, 50, 0.2, CrossBinary1(), 200, 10000, 3, 20, pr=False)
				s = g.work()
				sol_b.append(s.evaluate(f))

			cnt_d, cnt_b = 0, 0
			for s in sol_d:
				if s < EPS: cnt_d += 1
			for s in sol_b:
				if s < EPS: cnt_b += 1
			print "Double:", cnt_d, numpy.median(numpy.array(sol_d))
			print "Binary:", cnt_b, numpy.median(numpy.array(sol_b))
			print("----------------------------")

def zad4():
	best_p, best_m = 0, 0
	f = F6()
	for p in [30, 50, 100, 200]:
		sol = []
		for i in range(10):
			g = GeneticDouble(f, 2, -50, 150, 0.2, CrossDouble1(), p, 10000, 3, pr=False)
			s = g.work()
			sol.append(s.evaluate(f))
		m = numpy.median(numpy.array(sol))
		if best_p == 0 or m < best_m:
			best_p = p
			best_m = m
		print "P:", p, "| median:", m
	
	best_p_m, best_m = 0, 0
	for p in [0.1, 0.3, 0.5, 0.9]:
		sol = []
		for i in range(10):
			g = GeneticDouble(f, 2, -50, 150, p, CrossDouble1(), best_p, 10000, 3, pr=False)
			s = g.work()
			sol.append(s.evaluate(f))
		m = numpy.median(numpy.array(sol))
		print sol
		if best_p_m == 0 or m < best_m:
			best_p_m = p
			best_m = m
		print "p_m:", p, "| median:", m

	print "P:", best_p, "| p_m:", best_p_m

def zad5():
	f = F7()
	for t in [3, 5, 10, 20, 50]:
		sol = []
		for i in range(10):
			g = GeneticDouble(f, 2, -50, 150, 0.2, CrossDouble1(), 200, 10000, t, pr=False)
			s = g.work()
			sol.append(s.evaluate(f))
		m = numpy.median(numpy.array(sol))
		print "t_size:", t, "| median:", m

if argv[1] == "1" : zad1()
if argv[1] == "2" : zad2()
if argv[1] == "3" : zad3()
if argv[1] == "4" : zad4()
if argv[1] == "5" : zad5()