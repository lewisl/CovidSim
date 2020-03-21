#! /usr/bin/env python2.7
# -*- encoding: utf-8 -*-

import math

class Location(object):

    def __init__(self, x, y):
        """x and y are floats"""
        self.x = x
        self.y = y

    def move(self, deltaX, deltaY):
        """deltaX and deltaY are floats"""
        return Location(self.x + deltaX, self.y + deltaY)

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def distFrom(self, other):
        ox = other.x
        oy = other.y
        xDist = self.x - ox
        yDist = self.y - oy
        return (xDist**2 + yDist**2)**0.5

    def __str__(self):
        return '<' + str(self.x) + ', ' + str(self.y) + '>'


class Field(object):

    def __init__(self):
        self.drunks = {}

    def addDrunk(self, drunk, loc):
        if drunk in self.drunks:
            raise ValueError('Duplicate drunk')
        else:
            self.drunks[drunk] = loc

    def moveDrunk(self, drunk):
        if not drunk in self.drunks:
            raise ValueError('Drunk not in field')
        xDist, yDist = drunk.takeStep()
        currentLocation = self.drunks[drunk]
        #use move method of Location to get new location
        self.drunks[drunk] = currentLocation.move(xDist, yDist)

    def getLoc(self, drunk):
        if not drunk in self.drunks:
            raise ValueError('Drunk not in field')
        return self.drunks[drunk]

import random

class Drunk(object):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return 'This drunk is named ' + self.name

## the following are mine or downloaded
# class UsualDrunk(Drunk):
#     def takeStep(self):
#         stepChoices =\
#             [(0.0,1.0), (0.0,-1.0), (1.0, 0.0), (-1.0, 0.0)]
#         return random.choice(stepChoices)
#
#
# class ColdDrunk(Drunk):
#     def takeStep(self):
#         stepChoices =\
#             [(0.0,0.95), (0.0,-1.0), (1.0, 0.0), (-1.0, 0.0)]
#         return random.choice(stepChoices)
#
#
# class EDrunk(Drunk):
#     def takeStep(self):
#         deltaX = random.random()
#         if random.random() < 0.5:
#             deltaX = -deltaX
#         deltaY = random.random()
#         if random.random() < 0.5:
#             deltaY = -deltaY
#         return (deltaX, deltaY)


# the following are from the quiz
class UsualDrunk(Drunk):
    def takeStep(self):
        stepChoices =\
            [(0.0,1.0), (0.0,-1.0), (1.0, 0.0), (-1.0, 0.0)]
        return random.choice(stepChoices)

class ColdDrunk(Drunk):
    def takeStep(self):
        stepChoices =\
            [(0.0,0.9), (0.0,-1.03), (1.03, 0.0), (-1.03, 0.0)]
        return random.choice(stepChoices)

class EDrunk(Drunk):
    def takeStep(self):
        ang = 2 * math.pi * random.random()
        length = 0.5 + 0.5 * random.random()
        return (length * math.sin(ang), length * math.cos(ang))

class PhotoDrunk(Drunk):
    def takeStep(self):
        stepChoices =\
                    [(0.0, 0.5),(0.0, -0.5),
                     (1.5, 0.0),(-1.5, 0.0)]
        return random.choice(stepChoices)

class DDrunk(Drunk):
    def takeStep(self):
        stepChoices =\
                    [(0.85, 0.85), (-0.85, -0.85),
                     (-0.56, 0.56), (0.56, -0.56)]
        return random.choice(stepChoices)


def walk(f, d, numSteps):
    start = f.getLoc(d)
    for s in range(numSteps):
        f.moveDrunk(d)
    return(start.distFrom(f.getLoc(d)))

def walkVector(f, d, numSteps):
    start = f.getLoc(d)
    for s in range(numSteps):
        f.moveDrunk(d)
    return(f.getLoc(d).getX() - start.getX(),
           f.getLoc(d).getY() - start.getY())


def simWalks(numSteps, numTrials, dClass):
    homer = dClass('Homer')
    origin = Location(0, 0)
    distances = []
    for t in range(numTrials):
        f = Field()
        f.addDrunk(homer, origin)
        distances.append(walk(f, homer, numSteps))
    return distances

def simScatter(numSteps, numTrials, dClass):
    homer = dClass('Homer')
    origin = Location(0, 0)
    distancesX = []
    distancesY = []
    for t in range(numTrials):
        f = Field()
        f.addDrunk(homer, origin)
        dist = walkVector(f, homer, numSteps)
        distancesX.append(dist[0])
        distancesY.append(dist[1])
    return distancesX, distancesY

def drunkTest(numTrials = 20):
    random.seed(0)
    for numSteps in [10, 100, 1000, 10000]:
    # for numSteps in [0, 1]:
        distances = simWalks(numSteps, numTrials, EDrunk)
        print('Random walk of ' + str(numSteps) + ' steps')
        print(' Mean =', sum(distances)/len(distances))
        print(' Max =', max(distances), 'Min =', min(distances))

import pylab

def drunkTestP(numTrials = 50):
    stepsTaken = [10, 100, 1000, 10000]
    for dClass in (UsualDrunk, EDrunk, ColdDrunk):
        meanDistances = []
        for numSteps in stepsTaken:
            distances = simWalks(numSteps, numTrials, dClass)
            meanDistances.append(sum(distances)/len(distances))
        pylab.plot(stepsTaken, meanDistances, label = dClass.__name__)
        pylab.title('Mean Distance from Origin')
        pylab.xlabel('Steps Taken')
        pylab.ylabel('Steps from Origin')
        pylab.legend(loc = 'upper left')
        pylab.show()

def drunkScatter(dClass, numSteps = 100, numTrials = 50):
    distancesX,distancesY = simScatter(numSteps, numTrials, dClass)
    pylab.scatter(distancesX, distancesY, label = dClass.__name__)
    pylab.title('Distance from Origin')
    pylab.xlabel('X direction')
    pylab.ylabel('Y direction')
    pylab.legend(loc = 'upper left')
    pylab.show()

def drunkTestP1(numTrials = 50):
    stepsTaken = [10, 100, 1000, 10000]
    meanDistances = []
    squareRootOfSteps = []
    for numSteps in stepsTaken:
        distances = simWalks(numSteps, numTrials, EDrunk)
        meanDistances.append(sum(distances)/len(distances))
        squareRootOfSteps.append(numSteps**0.5)
    pylab.plot(stepsTaken, meanDistances, 'b-',
               label = 'Mean distance')
    pylab.plot(stepsTaken, squareRootOfSteps, 'g-.',
               label = 'Square root of steps')
    pylab.title('Mean Distance from Origin')
    pylab.xlabel('Steps Taken')
    pylab.ylabel('Steps from Origin')
    pylab.legend()
    pylab.show()