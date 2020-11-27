#Package requies numpy and matplotlib
#Install them with pip install numpy, and python -m pip install -U matplotlib

#Press c to toggle unit circle
#Press g to toggle gridlines
#Click on fundamental domain edge to highlight it and the paired side
#Use left/right arrow keys to move between consecutive sides
#Click on geodesic segment to highlight it and add an arrow
#Use up/down to move between consecutive segments
#Press a when a segment is highlighted to animate the drawing of the segment

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import path
from matplotlib import widgets
from matplotlib import collections
import matplotlib.animation as animation
from matplotlib.backend_bases import MouseButton
import time
import math
import sys
#plt.gcf().gca().add_artist(circle1)
#.gcf means get current figure, .gca means get current axis. Can also directly call plt.gca().

class animated_geodesic(object):
	def __init__(self, geod, op_man):#Init
		self.geod=geod
		self.colour=geod.colour
		self.nsides=geod.nsides
		self.op_man=op_man
		
		self.linewidth=4
		self.partsperside=20
		self.nframes=self.partsperside*len(self.geod.sides)
		self.isdone=0
		self.init_side(geod.curind) #Initialize the first side
		self.curpartofside=0 #Fix used for the first side, since the method is called twice.
		self.glinewidth=geod.linewidth
		
		self.finalind=self.curind #When we know we are done
		self.gsegs=[0]#0 is a placeholder to allow self.gsegs[-1]=... to work
		
		geod.coll.set_alpha(0)#Make invisible for now
		geod.removehighlight(True)
	
	def init_side(self, ind):#Initialize the constants for the current side
		self.curpartofside=1
		self.curind=ind #The current side
		self.curside=self.geod.sides[ind]
		self.theta1=self.geod.theta1s[ind]
		self.theta2=self.geod.theta2s[ind]
		self.dir=self.geod.directions[ind]
		if self.dir==1:#Forward arc
			self.theta=self.theta1
			self.isseg=False
		elif self.dir==-1:#Backwards arc
			self.theta=self.theta2
			self.isseg=False
		else:#We have a segment!
			self.theta=self.theta1
			self.isseg=True
		if self.isseg:
			self.angle=[(self.theta2[0]-self.theta1[0])/self.partsperside, (self.theta2[1]-self.theta1[1])/self.partsperside]#The dx, dy we use
		else:
			self.angle=(self.theta2-self.theta1)/self.partsperside
		self.append=True
	
	def animate(self, frames):
		if self.isdone:
			return []
		#First we check if we are moving onto a new side or continuing on the existing one.
		self.curpartofside+=1
		if self.curpartofside>self.partsperside:#Move on to new side
			ind=(self.curind+1)%self.nsides
			if ind==self.finalind:#Done!
				self.geod.coll.set_alpha(1)#Make visible again
				self.geod.ax.figure.canvas.draw()
				self.isdone=1
				return []
			self.init_side(ind) #Initialize the new side
		else:#Still on same side
			self.append=False
		
		side=self.curside #The side we are working on
		if self.isseg:
			self.theta=[self.theta[i]+self.angle[i] for i in range(2)]
			C=self.makeseg(self.theta1[0], self.theta1[1], self.theta[0], self.theta[1])
		else:
			if self.dir==1:#Normal, move forward
				self.theta+=self.angle
				C=self.makearc(self.geod.centres[self.curind][0], self.geod.centres[self.curind][1], self.geod.radii[self.curind], self.theta1, self.theta)
			else:#Backwards
				self.theta-=self.angle
				C=self.makearc(self.geod.centres[self.curind][0], self.geod.centres[self.curind][1], self.geod.radii[self.curind], self.theta, self.theta2)
		if self.append:#New side
			if len(self.gsegs)>0:
				self.gsegs[-1].set_linewidth(self.glinewidth)
			self.gsegs.append(self.geod.ax.add_patch(C))
		else:#Old side
			self.gsegs[-1]=self.geod.ax.add_patch(C)
		return self.gsegs

	def makearc(self, xcentre, ycentre, radius, startangle, endangle):#Make an arc
		return patches.Arc((xcentre, ycentre), 2*radius, 2*radius, angle=0, theta1=startangle, theta2=endangle, linewidth=self.linewidth, color=self.colour, zorder=1)

	def makeseg(self, startx, starty, endx, endy):#Make a segment
		return patches.Polygon(np.array([[startx, starty], [endx, endy]]), edgecolor=self.colour, facecolor='none', linewidth=self.linewidth, closed=False)

#Store funddom as [side pairing], then lines of xcentre, ycentre, radius, startangle, endangle	
class funddom(object):
	def __init__(self, ax, file, unitcirc):
		self.ax=ax #Storing axis
		self.linewidth=1
		self.colour='g'
		f=open(file, 'r') #Open it to read
		self.sidepairing=list(map(int, f.readline().split())) #Reading side pairing
		self.sides=[]
		for line in f:#Adding the sides of the fundamental domain in
			xcentre, ycentre, radius, iang1, iang2=map(float, line.split())#Getting data
			side=patches.Circle((xcentre, ycentre), radius)
			trans=side.get_transform() #To transform to the side
			arc=path.Path.arc(iang2, iang1) #Arc on unit circle
			self.sides.append(patches.PathPatch(trans.transform_path(arc))) #Adding the transformed arced side
		f.close()
		self.nsides=len(self.sides)#Number of sides
		self.coll = collections.PatchCollection(self.sides, edgecolor=self.colour, facecolor='none', linewidth=self.linewidth, zorder=1, picker=10)
		ax.add_collection(self.coll)
		self.hind=-1
		self.hoppind=-1

	def highlightpairing(self, ind):#Highlight paired edges
		if self.hind==ind:#Just return
			self.removehighlight(1)
			return
		self.removehighlight(0)
		otherside=self.sidepairing[ind]-1
		side=self.sides[ind]
		if ind==otherside:#Same side, make navy
			side.set_edgecolor("navy")
			side.set_facecolor("none")
			side.set_linewidth(4)
			self.ax.add_patch(side)
		else:#Different sides: make first side red, second side navy
			side.set_edgecolor("#ff0000")
			side.set_facecolor("none")
			side.set_linewidth(4)
			self.ax.add_patch(side)
			side=self.sides[otherside]
			side.set_edgecolor("navy")
			side.set_facecolor("none")
			side.set_linewidth(4)
			self.ax.add_patch(side)
		self.hind=ind
		self.hoppind=otherside
		self.ax.figure.canvas.draw()
	
	def removehighlight(self, update):#Remove the highlighting
		for i in [self.hind, self.hoppind]:
			try:
				self.sides[i].remove()
			except (ValueError, NotImplementedError):#If self.hind=self.hoppind or they =-1
				continue
		self.hind=-1
		self.hoppind=-1
		if update:
			self.ax.figure.canvas.draw()

#geodesic stores an entire geodesic, i.e. a set of arcs.
class geodesic(object):
	def __init__(self, ax, file, colour):#Init
		self.ax=ax #Storing axis
		self.colour=colour #Storing colour
		self.linewidth=1.3
		f=open(file, 'r') #Open it to read
		self.sides=[] #Stores the sides of the geodesic in order
		self.directions=[] #stores the directions of each segment
		self.theta1s=[]
		self.theta2s=[]
		self.radii=[]
		self.centres=[]
		for line in f:#Adding the sides of the fundamental domain in
			data=[float(dat) for dat in line.split()]
			if data[0]==0:#We have a circle arc. Data is stored as "0 xcentre ycentre radius startangle endangle direction"
				side=patches.Circle((data[1], data[2]), data[3])#Make the circle
				trans=side.get_transform() #The transformation
				arc=path.Path.arc(data[4], data[5])
				self.sides.append(patches.PathPatch(trans.transform_path(arc))) #Adding the transformed arced side
				self.directions.append(int(data[6]))
				self.theta1s.append(data[4])
				self.theta2s.append(data[5])
				self.radii.append(data[3])
				self.centres.append([data[1], data[2]])
			else:#We have a segment. Data is stored as "0 startx starty endx endy"
				start=[data[1], data[2]]
				end=[data[3], data[4]]
				self.sides.append(patches.Polygon(np.array([start, end]), closed=False))
				self.directions.append(int(0))
				self.theta1s.append(start)
				self.theta2s.append(end)
				self.radii.append(False)
				self.centres.append(False)
		f.close()

		#Add the sides to the collection
		self.nsides=len(self.sides)#Number of sides
		self.coll = collections.PatchCollection(self.sides, edgecolor=self.colour, facecolor='none', linewidth=self.linewidth, zorder=1, picker=8)
		ax.add_collection(self.coll)
			
		#Now constants relating to the highlighting of sides
		self.curind=-1
		self.arrow=0

	def makearrow(self, ind):
		side=self.sides[ind]
		verts=np.zeros((3, 2))#Stores the vertices
		r=self.radii[ind]
		if r!=False:
			theta1=self.theta1s[ind]
			theta2=self.theta2s[ind]
			if r>=0.5:
				totips=0.08/r
				halfbaselen=0.08
			else:
				totips=0.16
				halfbaselen=totips*r
			angle=(theta1+theta2)/2
			if theta1>theta2:#In case we loop past 0
				angle=angle+180
			angle=math.radians(angle)
			if self.directions[ind]==1:
				angle1=angle+totips
				angle2=angle-totips
			else:
				angle1=angle-totips
				angle2=angle+totips
			x, y=self.centres[ind]
			verts[0, 0]=x+math.cos(angle1)*r
			verts[0, 1]=y+math.sin(angle1)*r
			r2=r-halfbaselen
			verts[1, 0]=x+math.cos(angle2)*r2
			verts[1, 1]=y+math.sin(angle2)*r2
			r2=r+halfbaselen
			verts[2, 0]=x+math.cos(angle2)*r2
			verts[2, 1]=y+math.sin(angle2)*r2
			self.arrow=patches.Polygon(verts, color='black')
		else:
			start=self.theta1s[ind]
			end=self.theta2s[ind]
			dx=end[0]-start[0]
			dy=end[1]-start[1]
			hypot=dx**2+dy**2
			if hypot<=1:
				totips=0.16*hypot
				halfbaselen=0.08*hypot
			else:
				totips=0.16
				halfbaselen=0.08
			dx=dx/hypot
			dy=dy/hypot
			middle=[(start[i]+end[i])/2 for i in range(2)]
			verts[0,]=[middle[0]+dx*totips, middle[1]+dy*totips]
			dx1=dx*halfbaselen
			dy1=dy*halfbaselen
			verts[1,]=[middle[0]-dy1, middle[1]+dx1]
			verts[2,]=[middle[0]+dy1, middle[1]-dx1]
			self.arrow=patches.Polygon(verts, color='black')

	#Highlight a segment
	def highlight(self, ind):
		ind=ind%self.nsides
		if self.curind==ind:#Just return
			self.removehighlight(1)
			return 0
		self.removehighlight(0)
		i=ind-1
		c=['crimson', 'orange']
		for j in range(2):
			i=(i+1)%self.nsides
			self.sides[i].set_edgecolor(c[j])
			self.sides[i].set_facecolor("none")
			self.sides[i].set_linewidth(3)
			self.ax.add_patch(self.sides[i])
		self.curind=ind
		self.makearrow(ind)
		self.ax.add_patch(self.arrow)
		self.ax.figure.canvas.draw()
		return 1
	
	def removehighlight(self, update):#Remove the highlighting
		for i in [self.curind, (self.curind+1)%self.nsides]:
			try:
				self.sides[i].remove()
			except (ValueError, NotImplementedError):
				continue
		self.curind=-1
		try:
			self.arrow.remove()
		except (ValueError, AttributeError):
			pass
		if update:
			self.ax.figure.canvas.draw()

#Handles the user input
class operations_manager(object):
	def __init__(self, ax, fdom, geodesics, unitcirc):#Storing them for now
		self.ax=ax
		self.fdom=fdom
		self.geodesics=geodesics
		self.unitcirc=unitcirc
		self.geocoll=[geod.coll for geod in geodesics]#The collections of arcs/segments in the geodesics
		self.hgeod=-1
	
	def keypress(self, event):
		if self.fdom!=0:
			ind=self.fdom.hind
			if ind!=-1:#We want to move left or right
				if event.key=='right':
					ind=(ind+1)%self.fdom.nsides
					self.fdom.highlightpairing(ind)
				elif event.key=='left':
					ind=(ind-1)%self.fdom.nsides
					self.fdom.highlightpairing(ind)
		if self.hgeod!=-1:
			if event.key=='up':
				self.geodesics[self.hgeod].highlight(self.geodesics[self.hgeod].curind+1)
			elif event.key=='down':
				self.geodesics[self.hgeod].highlight(self.geodesics[self.hgeod].curind-1)
			elif event.key=='a':
				anim_geod=animated_geodesic(self.geodesics[self.hgeod], self)
				anim=animation.FuncAnimation(self.ax.figure, anim_geod.animate, frames=anim_geod.nframes, interval=25, blit=True, repeat=False)
				self.hgeod=-1
				plt.show()
			#elif event.key=='A':#For saving the animated geodesic
				#anim_geod=animated_geodesic(self.geodesics[self.hgeod], self)
				#anim=animation.FuncAnimation(self.ax.figure, anim_geod.animate, frames=anim_geod.nframes, interval=25, blit=True, repeat=False)
				#anim.save('testinganim.gif')
				#plt.show()
		if event.key=='c':
			self.unitcirc.set_alpha(1-self.unitcirc.get_alpha())
			self.ax.figure.canvas.draw()
	
	def onclick(self, event):
		if event.mouseevent.button!=MouseButton.LEFT:
			return
		if event.artist==self.fdom.coll:#Clicked in fundamental domain
			self.fdom.highlightpairing(event.ind[0])
		else:
			try:
				ind=self.geocoll.index(event.artist)
			except ValueError:
				return
			for i in event.ind:
				if ind!=self.hgeod:#Removing highlighting on other geodesics
					self.geodesics[self.hgeod].removehighlight(1)
				if self.geodesics[ind].highlight(i):#Storing the highlighted geodesic
					self.hgeod=ind
				else:#We were removing the highlighting
					self.hgeod=-1
				return

#The main body
if __name__=='__main__':
	#Initializing the figure
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	plt.axis('square')
	fig.tight_layout()
	
	#Axes setting
	xymin=-1.03
	xymax=1.030
	xytickd=0.25
	skip=2
	major_ticks = np.arange(-1, 1.1, skip*xytickd)
	minor_ticks = np.arange(-1, 1.1, xytickd)
	ax.set_xlim(xymin, xymax)
	ax.set_xticks(major_ticks)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_ylim(xymin, xymax)
	ax.set_yticks(major_ticks)
	ax.set_yticks(minor_ticks, minor=True)
	
	#The unit circle
	unitcirc=patches.Circle((0, 0), 1, linewidth=1, edgecolor="b", facecolor="none", zorder=1, picker=1, alpha=1)
	ax.add_patch(unitcirc)
	
	#Geodesic arcs
	geodesics=[]
	arccolours=['aqua', 'chartreuse', 'fuchsia', 'blue', 'maroon', 'darkblue', 'teal', 'y']
	counter=-1
	isfdom=False
	for c in sys.argv[1:]:
		if c[0]=="f": #Fundamental domain
			fdom=funddom(ax,'fdoms/'+c+'.dat', unitcirc)
			isfdom=True
		else:
			counter=(counter+1)%len(arccolours)
			geodesics.append(geodesic(ax,'fdoms/'+c+'.dat', arccolours[counter]))
	if isfdom:
		op_man=operations_manager(ax, fdom, geodesics, unitcirc)
	else:
		op_man=operations_manager(ax, 0, geodesics, unitcirc)
	
	matplotlib.rcParams['keymap.back'].remove('left')
	matplotlib.rcParams['keymap.forward'].remove('right')
	fig.canvas.mpl_connect('pick_event', op_man.onclick)
	fig.canvas.mpl_connect('key_press_event', op_man.keypress)
	plt.show()