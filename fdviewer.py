#Package requies numpy and matplotlib
#Install them with pip install numpy, and python -m pip install -U matplotlib

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import widgets
import math
import sys
#from matplotlib.widgets import Button
#from matplotlib.widgets import RadioButtons

#plt.gcf().gca().add_artist(circle1)
#.gcf means get current figure, .gca means get current axis. Can also directly call plt.gca().

#arc stores a single arc to plot.
class arc(object):
	def __init__(self, ax, xcentre, ycentre, radius, startangle, endangle, direction, colour):
		self.seg=ax.add_patch(patches.Arc((xcentre, ycentre), 2*radius, 2*radius, angle=0, theta1=startangle, theta2=endangle, linewidth=1.3, color=colour, zorder=1)) #The arc
		self.dir=direction #The direction
	
	def animate(self):#TO DO
		print(1)
	
	def change_colour(self, c):#TO DO
		print("colour changed")

class segment(object):
	def __init__(self, ax, startx, starty, endx, endy, colour):
		self.ax=ax #Store axis
		self.xpath=[startx, endx]
		self.ypath=[starty, endy]
		self.seg=ax.plot(self.xpath, self.ypath, color=colour)

#geodesic stores an entire geodesic, i.e. a set of arcs.
class geodesic(object):
	def __init__(self, ax, file, colour):
		self.ax=ax #Storing axis
		f=open(file, 'r') #Open it to read
		self.sides=[] #Stores the sides of the geodesic in order
		for line in f:#Adding the sides of the fundamental domain in
			data=[float(dat) for dat in line.split()]
			if data[0]==0:#We have a circle arc. Data is stored as "0 xcentre ycentre radius startangle endangle direction"
				self.sides.append(arc(self.ax, data[1], data[2], data[3], data[4], data[5], data[6], colour))
			else:#We have a segment. Data is stored as "0 startx starty endx endy"
				self.sides.append(segment(self.ax, data[1], data[2], data[3], data[4], colour))
		f.close()
	
#Store funddom as [side pairing], then lines of xcentre, ycentre, radius, startangle, endangle	
class funddom(object):
	def __init__(self, ax, file):
		self.ax=ax #Storing axis
		f=open(file, 'r') #Open it to read
		self.sidepairing=list(map(float, f.readline().split())) #Reading side pairing
		self.sides=[]
		self.vangles=list(map(float, f.readline().split()))#Reading vertex angles
		for line in f:#Adding the sides of the fundamental domain in
			xcentre, ycentre, radius, startangle, endangle=map(float, line.split())#Getting data
			self.sides.append(arc(self.ax, xcentre, ycentre, radius, startangle, endangle, 0, "g"))#No direction!
			#Let's add in the middle angle!
		f.close()

class operations_manager(object):
	def __init__(self, ax, fdom, geodesics):#Storing them for now
		self.ax=ax
		self.fdom=fdom
		self.geodesics=geodesics
	
#tracker for clicking on axis. Will highlight the corresponding side pairing
#Buttons to animate each arc.	
	
if __name__=='__main__':
	#Initializing the figure
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	plt.axis('scaled')
	
	#Axes setting
	xymin=-1.2
	xymax=1.201
	xytickd=0.2
	skip=3
	major_ticks = np.arange(xymin, xymax, skip*xytickd)
	minor_ticks = np.arange(xymin, xymax, xytickd)
	ax.set_xlim(xymin, xymax)
	ax.set_xticks(major_ticks)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_ylim(xymin, xymax)
	ax.set_yticks(major_ticks)
	ax.set_yticks(minor_ticks, minor=True)
	ax.grid(which='both', markeredgewidth=0.1, zorder=0)
	
	#The unit circle
	arc(ax, 0, 0, 1, 0, 360, 0, 'b')
	
	#Geodesic arcs
	geodesics=[]
	arccolours=['r', 'm', 'y', 'c']
	counter=-1
	isfdom=False
	for c in sys.argv[1:]:
		if c[0]=="f": #Fundamental domain
			fdom=funddom(ax,'fdoms/'+c+'.dat')
			isfdom=True
		else:
			counter=(counter+1)%len(arccolours)
			geodesics.append(geodesic(ax,'fdoms/'+c+'.dat', arccolours[counter]))
	if isfdom:
		op_man=operations_manager(ax, fdom, geodesics)
	else:
		op_man=operations_manager(ax, 0, geodesics)
		
	#xrange=widgets.Slider(ax,'width', 0.1, 10, valinit=1.2, valfmt='%1.2f', closedmin=True, closedmax=True, slidermin=None, slidermax=None, dragging=True, valstep=None, orientation='horizontal')	
	#xrange.on_changed(update(ax,xrange.val))
	
	plt.show()