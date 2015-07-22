#! /usr/bin/env/python

"""Softward tools to simulate making etchings and printing them on paper.
claw
4 March 2012
"""

import pylab as p
import numpy as n
import matplotlib.colors as colors

def point_inside_polygon(x,y,poly):
    """Needed to convert polygon etchings into bitmapped masks for printing on paper.
    """

    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

class etching():
    """Defines etching object. This is basically a single block (of say linoleum) with axes from 0 to 1.
    Each etch removes material and creates a boolean mask marked False in etched area.
    This can be initialized with a mask from previous etching. Etch method then returns only the newly etched area.
    """

    def __init__(self, mask=[0]):
        """Get new etching surface
        """

        self.size = 300
        p.figure(1)
        p.clf()
        p.axis([0,1,0,1])
        if len(mask) == 1:
            self.mask = n.ones((self.size,self.size), dtype='bool')
        else:
            self.mask = mask
        print 'Surface ready.'

    def etch(self, refmask=[0]):
        """Cut away at surface. Can only click polygons for now.
        Optionally input reference mask to guide etching.
        """

        size = self.size
        mask = n.ones((self.size,self.size), dtype='bool')

        print 'Click for points of region to etch (right click to exit).'
        p.figure(1)
        if len(refmask) != 1:
            p.imshow(n.transpose(-refmask), aspect='auto', origin='lower', interpolation='nearest', cmap=p.cm.Greys, extent=(0,1,0,1), vmax=0.5)
        else:
            p.imshow(n.transpose(-self.mask), aspect='auto', origin='lower', interpolation='nearest', cmap=p.cm.Greys, extent=(0,1,0,1), vmax=0.5)
        xy = p.ginput(n=0, timeout=3000)
        xy = n.array(xy)
        
        print 'Calculating remaining region...'
        p.figure(1)
        p.axis([0,1,0,1])
        p.fill(xy[:,0],xy[:,1],'k')

        for i in range(size):
            for j in range(size):
                mask[i,j] = not(point_inside_polygon(float(i)/size, float(j)/size, xy))
        self.mask = self.mask * mask


class paper():
    """Piece of paper to print etchings on. Extends from 0 to 1.
    Each press has a color that covers previous color in etching mask area that is True.
    """

    def __init__(self):
        """Get fresh piece of paper.
        """

        self.size = 300
        p.figure(2)
        p.clf()
        p.axis([0,1,0,1])
        self.clist = 'w'
        self.mask = n.zeros((self.size,self.size), dtype='int')
        print 'Paper ready.'

    def press(self, etching, color='k'):
        """Press etch onto paper. Optionally can define a color of ink.
        """

        p.figure(2)
        p.axis([0,1,0,1])

        self.clist = self.clist + color
        maskval = len(self.clist)-1  # set mask to this for this color
        cmap = colors.ListedColormap(self.clist)
        self.mask = n.where(etching.mask, maskval*n.ones((self.size, self.size), dtype='int'), self.mask)

        p.imshow(n.transpose(self.mask), aspect='auto', origin='lower', interpolation='nearest', cmap=cmap, extent=(0,1,0,1))
        print 'Pressed etching onto paper.'
