#!/usr/bin/env python

"""shors.py: Shor's algorithm for quantum integer factorization"""

import math
import random
import argparse
import threading
import time

__author__ = "Todd Wildey"
__copyright__ = "Copyright 2013"
__credits__ = ["Todd Wildey"]

__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Todd Wildey"
__email__ = "toddwildey@gmail.com"
__status__ = "Prototype"

def printNone(str):
	pass

def printVerbose(str):
	print(str)

printInfo = printNone

####################################################################################################
#                                                                                                   
#                                        Quantum Components                                         
#                                                                                                   
####################################################################################################

class Mapping:
	def __init__(self, state, amplitude):
		self.state = state
		self.amplitude = amplitude


class QuantumState:
	def __init__(self, amplitude, register):
		self.amplitude = amplitude
		self.register = register
		self.entangled = {}

	def entangle(self, fromState, amplitude):
		register = fromState.register
		entanglement = Mapping(fromState, amplitude)
		try:
			self.entangled[register].append(entanglement)
		except KeyError:
			self.entangled[register] = [entanglement]

	def entangles(self, register = None):
		entangles = 0
		if register is None:
			for states in self.entangled.values():
				entangles += len(states)
		else:
			entangles = len(self.entangled[register])

		return entangles


class QubitRegister:
	def __init__(self, numBits):
		self.numBits = numBits
		self.numStates = 1 << numBits
		self.entangled = []
		
		print(str("+1\t" + str(self.numStates)))
		self.states = [QuantumState(complex(0.0), self) for x in range(self.numStates)]
		print(str("+2"))
		self.states[0].amplitude = complex(1.0)
		print(str("+3"))

	def propagate(self, fromRegister = None):
		if fromRegister is not None:
			for state in self.states:
				amplitude = complex(0.0)

				try:
					entangles = state.entangled[fromRegister]
					for entangle in entangles:
						amplitude += entangle.state.amplitude * entangle.amplitude

					state.amplitude = amplitude
				except KeyError:
					state.amplitude = amplitude

		for register in self.entangled:
			if register is fromRegister:
				continue

			register.propagate(self)

	# Map will convert any mapping to a unitary tensor given each element v
	# returned by the mapping has the property v * v.conjugate() = 1
	#
	def map(self, toRegister, mapping, propagate = True):
		self.entangled.append(toRegister)
		toRegister.entangled.append(self)

		# Create the covariant/contravariant representations
		mapTensorX = {}
		mapTensorY = {}
		for x in range(self.numStates):
			mapTensorX[x] = {}
			codomain = mapping(x)
			for element in codomain:
				y = element.state
				mapTensorX[x][y] = element

				try:
					mapTensorY[y][x] = element
				except KeyError:
					mapTensorY[y] = { x: element }

		# Normalize the mapping:
		def normalize(tensor, p = False):
			lSqrt = math.sqrt
			for vectors in tensor.values():
				sumProb = 0.0
				for element in vectors.values():
					amplitude = element.amplitude
					sumProb += (amplitude * amplitude.conjugate()).real

				normalized = lSqrt(sumProb)
				for element in vectors.values():
					element.amplitude = element.amplitude / normalized

		normalize(mapTensorX)
		normalize(mapTensorY, True)

		# Entangle the registers
		for x, yStates in mapTensorX.items():
			for y, element in yStates.items():
				amplitude = element.amplitude
				toState = toRegister.states[y]
				fromState = self.states[x]
				toState.entangle(fromState, amplitude)
				fromState.entangle(toState, amplitude.conjugate())

		if propagate:
			toRegister.propagate(self)

	def measure(self):
		measure = random.random()
		sumProb = 0.0

		# Pick a state
		finalX = None
		finalState = None
		for x, state in enumerate(self.states):
			amplitude = state.amplitude
			sumProb += (amplitude * amplitude.conjugate()).real

			if sumProb > measure:
				finalState = state
				finalX = x
				break

		# If state was found, update the system
		if finalState is not None:
			for state in self.states:
				state.amplitude = complex(0.0)

			finalState.amplitude = complex(1.0)
			self.propagate()

		return finalX

	def entangles(self, register = None):
		entangles = 0
		for state in self.states:
			entangles += state.entangles(None)

		return entangles

	def amplitudes(self):
		amplitudes = []
		for state in self.states:
			amplitudes.append(state.amplitude)

		return amplitudes

def printEntangles(register):
	printInfo("Entagles: " + str(register.entangles()))

def printAmplitudes(register):
	amplitudes = register.amplitudes()
	for x, amplitude in enumerate(amplitudes):
		printInfo('State #' + str(x) + '\'s amplitude: ' + str(amplitude))

def hadamard(x, Q):
	codomain = []
	for y in range(Q):
		amplitude = complex(pow(-1.0, bitCount(x & y) & 1))
		codomain.append(Mapping(y, amplitude))

	return  codomain

# Quantum Modular Exponentiation
def qModExp(a, exp, mod):
	state = modExp(a, exp, mod)
	amplitude = complex(1.0)
	return [Mapping(state, amplitude)]

# Quantum Fourier Transform
def qft(x, Q):
	fQ = float(Q)
	k = -2.0 * math.pi
	codomain = []

	for y in range(Q):
		theta = (k * float((x * y) % Q)) / fQ
		amplitude = complex(math.cos(theta), math.sin(theta))
		codomain.append(Mapping(y, amplitude))

	return codomain


def pow_mod(a, b, n):
	temp = 1
	while b > 0:
		temp = (temp * a) % n
		b -= 1
	return temp
	
Period = 0
d = 0 
def findPeriod_job(a, N, r_start, r_end):
	global d
	global Period
	r = r_start
	
	temp =  (a ** r) % N # pow_mod(a, r-1, N)
	print("findPeriod_job r_start\t" + str(r_start) + "\ttemp = " + str(temp))
	while(r < r_end and Period == 0):
		#print(str("findPeriod_job temp\t") + str(temp))
		temp = (temp * a) % N
		r += 1
			
		if temp == a:
			r -= 1
			d = r
			Period = 1
			print(str("found\t" + "r = " + str(r)))
			break
	

def findPeriod(a, N, pc = 2):

	if pc == 1:
		r = 1
		temp = a
		while(True):
			temp *= a 
			r += 1
			
			if(temp > N):
				temp %= N
				
			if temp == a:
				r -= 1
				break
				
	else:
		if pc == 2:
			max_treads = 14
			r_start = 2
			size = math.floor(a/max_treads)
			global d
			global Period
			while(max_treads > 0):
				t = threading.Thread(target=findPeriod_job, args=(a, N, r_start, r_start + size))
				t.start()
				r_start = r_start + size
				# print("r = \t" + str(r))
				max_treads -= 1 
			while Period == 0:
				pass
			Period = 0
			print("d = " + str(d))
			return d

		else:
			if pc == 3:
				r = math.floor(math.sqrt(N))
				temp = a ** r
				while(True):
					# print(str(temp))
					if(temp > N):
						temp %= N
						
					if temp == a:
						r -= 1
						break
					r -= 1
					temp = a ** r
			else:
				nNumBits = N.bit_length()
				inputNumBits = (2 * nNumBits) - 1
				inputNumBits += 1 if ((1 << inputNumBits) < (N * N)) else 0
				Q = 1 << inputNumBits

				printInfo("Finding the period...")
				printInfo("Q = " + str(Q) + "\ta = " + str(a))
				printInfo("1")
				inputRegister = QubitRegister(inputNumBits)
				printInfo("1.1")
				hmdInputRegister = QubitRegister(inputNumBits)
				printInfo("1.2")
				qftInputRegister = QubitRegister(inputNumBits)
				printInfo("1.3")
				outputRegister = QubitRegister(inputNumBits)
				printInfo("2")
				printInfo("Registers generated")
				printInfo("Performing Hadamard on input register")

				inputRegister.map(hmdInputRegister, lambda x: hadamard(x, Q), False)
				# inputRegister.hadamard(False)
				printInfo("3")
				printInfo("Hadamard complete")
				printInfo("Mapping input register to output register, where f(x) is a^x mod N")

				hmdInputRegister.map(outputRegister, lambda x: qModExp(a, x, N), False)
				printInfo("4")
				printInfo("Modular exponentiation complete")
				printInfo("Performing quantum Fourier transform on output register")

				hmdInputRegister.map(qftInputRegister, lambda x: qft(x, Q), False)
				inputRegister.propagate()
				printInfo("5")
				printInfo("Quantum Fourier transform complete")
				printInfo("Performing a measurement on the output register")

				y = outputRegister.measure()

				printInfo("Output register measured\ty = " + str(y))

				# Interesting to watch - simply uncomment
				# printAmplitudes(inputRegister)
				# printAmplitudes(qftInputRegister)
				# printAmplitudes(outputRegister)
				# printEntangles(inputRegister)

				printInfo("Performing a measurement on the periodicity register")
				printInfo("6")
				x = qftInputRegister.measure()

				printInfo("QFT register measured\tx = " + str(x))

				if x is None:
					return None

				printInfo("Finding the period via continued fractions")
				printInfo("7")
				r = cf(x, Q, N)
				printInfo("8")
				printInfo("Candidate period\tr = " + str(r))

	return r

####################################################################################################
#                                                                                                   
#                                       Classical Components                                        
#                                                                                                   
####################################################################################################

BIT_LIMIT = 4096

def bitCount(x):
	sumBits = 0
	while x > 0:
		sumBits += x & 1
		x >>= 1

	return sumBits

# Greatest Common Divisor
def gcd(a, b):
	while b != 0:
		tA = a % b
		a = b
		b = tA

	return a

# Extended Euclidean
def extendedGCD(a, b):
	fractions = []
	while b != 0:
		fractions.append(a // b)
		tA = a % b
		a = b
		b = tA

	return fractions

# Continued Fractions
def cf(y, Q, N):
	fractions = extendedGCD(y, Q)
	depth = 2

	def partial(fractions, depth):
		c = 0
		r = 1

		for i in reversed(range(depth)):
			tR = fractions[i] * r + c
			c = r
			r = tR

		return c

	r = 0
	for d in range(depth, len(fractions) + 1):
		tR = partial(fractions, d)
		if tR == r or tR >= N:
			return r

		r = tR

	return r

# Modular Exponentiation
def modExp(a, exp, mod):
	fx = 1
	while exp > 0:
		if (exp & 1) == 1:
			fx = fx * a % mod
		a = (a * a) % mod
		exp = exp >> 1

	return fx

def pick(N):
	a = math.floor((random.random() * (N - 1)) + 0.5)
	return a

def checkCandidates(a, r, N, neighborhood):
	if r is None:
		return None

	# Check multiples
	for k in range(1, neighborhood + 2):
		tR = k * r
		if modExp(a, a, N) == modExp(a, a + tR, N):
			return tR

	# Check lower neighborhood
	for tR in range(r - neighborhood, r):
		if modExp(a, a, N) == modExp(a, a + tR, N):
			return tR

	# Check upper neigborhood
	for tR in range(r + 1, r + neighborhood + 1):
		if modExp(a, a, N) == modExp(a, a + tR, N):
			return tR

	return None

def shors(N, attempts = 1, neighborhood = 0.0, numPeriods = 1):
	if(N.bit_length() > BIT_LIMIT or N < 3):
		return False

	periods = []
	neighborhood = math.floor(N * neighborhood) + 1

	printInfo("N = " + str(N))
	printInfo("Neighborhood = " + str(neighborhood))
	printInfo("Number of periods = " + str(numPeriods))
	
	for attempt in range(attempts):
		

		a = pick(100)
		while a < 2:
			a = pick(100)
		
		d = gcd(a, N)
		if d != 1:
			printInfo("Found factors classically, re-attempt" + "\td = " + str((d)))
			continue

		printInfo("\nAttempt #" + str(attempt) + "\ta = " + str(a))
		r = findPeriod(a, N)
		
		printInfo("Checking candidate period, nearby values, and multiples")

		#r = checkCandidates(a, r, N, neighborhood)

		#if r is None:
		#	printInfo("Period was not found, re-attempt")
		#	continue

		if (r % 2) > 0:
			printInfo("Period was odd, re-attempt")
			continue

		d = modExp(a, (r // 2), N)
		if r == 0 or d == (N - 1):
			printInfo("Period was trivial, re-attempt\t" + str(r))
			continue
		
		printInfo("Period found\tr = " + str(r) + "\ta = " + str(a))

		periods.append(r)
		if(len(periods) < numPeriods):
			continue

		printInfo("\nFinding least common multiple of all periods")

		r = 1
		for period in periods:
			d = gcd(period, r)
			r = (r * period) // d

		b = modExp(a, (r // 2), N)
		f1 = gcd(N, b + 1)
		f2 = gcd(N, b - 1)

		return [f1, f2]

	return None

####################################################################################################
#                                                                                                   
#                                    Command-line functionality                                     
#                                                                                                   
####################################################################################################

def parseArgs():
	parser = argparse.ArgumentParser(description='Simulate Shor\'s algorithm for N.')
	parser.add_argument('-a', '--attempts', type=int, default=20, help='Number of quantum attemtps to perform')
	parser.add_argument('-n', '--neighborhood', type=float, default=0.01, help='Neighborhood size for checking candidates (as percentage of N)')
	parser.add_argument('-p', '--periods', type=int, default=2, help='Number of periods to get before determining least common multiple')
	parser.add_argument('-v', '--verbose', type=bool, default=True, help='Verbose')
	parser.add_argument('N', type=int, help='The integer to factor')
	return parser.parse_args()



# Example to n in 64 bits  4611686181636145867 = 2147483659*2147483713
# Example to n in 256 bits 61946200040551647909307739416522924093466643837971988782732316819163396610453 = 289798385584917605485498001977825965287 * 213756194381559050213318332733439210019

# Example to n in 512 bits 7481941016577287348099705569618881449446613674024887958060868042015962068162606112448282331657456711392343427317702862486293906907971066705762469196914429
    
# Example to n in 1024 bits 92562033715812759124114741736044283523381289941048932825510555468820581347335614003659313937583288939531366470214486407630811509181413157179309088380371985857285673607695187430424172796543432626204701421657401048700926260108171681237842250547938815679826497511311277132455416142287736205330601381198436096309


def main():
	args = parseArgs()

	global printInfo
	if args.verbose:
		printInfo = printVerbose
	else:
		printInfo = printNone

	factors = shors(args.N, args.attempts, args.neighborhood, args.periods)
	if factors is not None:
		print("Factors:\t" + str(factors[0]) + ", " + str(factors[1]))

if __name__ == "__main__":
	main()