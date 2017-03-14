##This script reads the result of save save_sift_correspondenc.ijm
## and fits a plane or a quadratic to the saved points
## Will reconstruct slice and full stack.

from PIL import Image
import os, glob
import numpy as np
import os.path
import sys
import math
import itertools
import time
import re
import functions as DIL
import numpy as np
import scipy.linalg
import scipy.optimize
from scipy import stats
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import pickle

def rotate_stack(astack,RotM):
	RotMinv = np.linalg.inv(RotM)
	xdim, ydim, zdim =np.shape(astack)
	icent,jcent,kcent = ((xdim+1)/2, (ydim+1)/2, (zdim+1)/2)
	max_dim = max([xdim,ydim,zdim])
	rotated = np.zeros((max_dim,max_dim,max_dim))
	for kstar in range(max_dim):
		for istar in range(max_dim):
			for jstar in range(max_dim):
				[i,j,k] = np.dot(RotMinv,np.array([istar-icent,jstar-jcent,kstar-kcent]))+np.array([icent,jcent,kcent])
				i=int(round(i))
				j=int(round(j))
				k=int(round(k))
				if i<xdim and i>=0 and j<ydim and j>=0 and k<zdim and k>=0: #make sure we are not out of bounds 
					rotated[istar,jstar,kstar] = astack[i,j,k] #round is nearest neighbour
				else:
					pass
	return rotated
	
			
##LOAD SIFT FIJI DATA
def load_fiji_sift(direct):
	#direct is the folder holding the txt files from fiji save_sift_correspondence.ijm
	#returns an array of points
	save_points = []
	
	for txt in sorted(glob.glob("{0}*.txt".format(direct))):
		numbers = re.findall('\d+',txt) #extract numbers from a string
		slice_no = int(numbers[-1]) #convert to int 					#CAREFULHERE we use last number in string as slice number
		info = np.loadtxt(txt, skiprows= 1)
		#print 'slice no is {0}'.format(slice_no)
		if len(info) ==5: # skip empty
			pass
		else:
			no_points = np.shape(info)[0]
			for i in range(no_points):
				point_row = info[i]
				point = [int(round(point_row[6])),int(round(point_row[5])),slice_no] #sticking to image axis convention
				save_points.append(point)
	
					
	data = np.array(save_points)
	return data

def plot_points(points,STACK,save_dir='./point_stack/'):
	#takes output from load_fiji_sift
	#and makes array of the points
	stackX,stackY,stackZ =np.shape(STACK)
	point_stack = np.zeros((stackX,stackY,stackZ))
	for point in points:
		point_stack[point[0],point[1],point[2]] = 255
	DIL.save_as_im_sequence(point_stack,save_dir)


def plane_fit(points,STACK,order=1,plot=False,slice_map=True,save_image=True, save_name='',save_stack=False,
		save_stack_dir=''):
	#fits plane to points if order = 1
	#returns C a vecor such that in the oprder 1 case
	#Z = C[0]*X +C[1]*Y + C[2]
	#in order 2 case
	#Z = C[0] + C[1]*X +C[2]*Y +C[3]*X*Y +C[4]*X**2 +C[5]Y**2
	#STACK is 3numpyarray
	#WHEN PLOTTING REDUCE THE SIZE OF X and Y !!!
	stackX,stackY,stackZ = np.shape(STACK)
	X,Y = np.meshgrid(np.arange(0, stackY,1), np.arange(0,stackX,1))
	XX = X.flatten()
	YY = Y.flatten()
	slice = np.zeros((stackX,stackY)) #to write too
	if save_stack == True:
		aligned_stack = np.empty(np.shape(STACK))
		aligned_stack[:] = -1
	if order == 1:
		# best-fit linear plane
		A = np.c_[points[:,0], points[:,1], np.ones(points.shape[0])]
    		C,_,_,_ = scipy.linalg.lstsq(A, points[:,2])    # coefficients
    
   		# evaluate it on grid
		Z = C[0]*X + C[1]*Y + C[2]
		if slice_map == True:
			slice_image = np.empty((stackX,stackY))
			slice_image[:] = -1
		#for save image:
		if save_image ==True:
			for i in range(stackX):
				for j in range(stackY):
					z= int(round(i*C[0] + j*C[1] + C[2]))
					try:
						slice[i,j] = STACK[i,j,z]
						slice_image[i,j] = z
						
					except IndexError:
						pass
		if save_stack == True:
			for x in range(stackX):
				for y in range(stackY):
					for z in range(stackZ):
						new_z= int(round(x*C[0] + y*C[1] + z))
						try:
							aligned_stack[x,y,z] = STACK[x,y,new_z]
						except IndexError:
							pass
					

		
	elif order == 2:
    		# best-fit quadratic curve
		A = np.c_[np.ones(points.shape[0]), points[:,:2], np.prod(points[:,:2], axis=1), points[:,:2]**2]
		C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
  		#evaluate	
		Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
		if slice_map == True:
			slice_image = np.empty((stackX,stackY))
			slice_image[:] = -1
		#for save
		if save_image == True:
			for i in range(stackX):
				for j in range(stackY):
					z= int(round(C[0] + i*C[1] + j*C[2]+j*i*C[3]+C[4]*i**2+C[5]*j**2))
					try:
						slice[i,j] = STACK[i,j,z]
						if slice_map == True:
							slice_image[i,j] = z
					except IndexError:
						pass
		if save_stack == True:
			for x in range(stackX):
				for y in range(stackY):
					for z in range(stackZ):
						new_z = int(round(z+ x*C[1] + y*C[2]+y*x*C[3]+C[4]*x**2+C[5]*y**2))
						try:
							aligned_stack[x,y,z] = STACK[x,y,new_z]

						except IndexError:
							pass



	if save_image == True:
		out = Image.fromarray(slice) #SOLUTION
		out.save('{0}_order{1}.tif'.format(save_name,order))
		out.close()
		if slice_map == True:
			slice_image_img = Image.fromarray(slice_image)
			slice_image_img.save('{0}_order{1}_slice_map.tif'.format(save_name,order))
			slice_image_img.close()

	if save_stack == True:
		DIL.save_as_im_sequence(aligned_stack,save_stack_dir)
	if plot==True:
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
		ax.scatter(points[:,0], points[:,1], points[:,2], c='r', s=50)
		plt.xlabel('X')
		plt.ylabel('Y')
		ax.set_zlabel('Z')
		ax.axis('equal')
		ax.axis('tight')
		plt.show()
	return C





def filter_out(tuples,param):
	tuplesx = [x[0] for x in tuples]
	tuplesy = [x[1] for x in tuples]
	xc = np.mean(tuplesx)
	yc = np.mean(tuplesy)
	centre = (xc,yc)
	filtered_points = []
	for tup in tuples:
		norm = ((xc-tup[0])**2+(yc-tup[1])**2)**0.5
		if norm < param:
			filtered_points.append(tup)

	return filtered_points

def points_inside(poly):
	#returns a list of points inside a polygon defined
	#by a list of vertices poly.
	tick = time.time()
	poly = np.array(poly)
	polypath = mplPath.Path(poly)
	polyx = poly[:,0]
	polyy= poly[:,1]
	xmax = polyx.max()
	xmin = polyx.min()
	ymax = polyy.max()
	ymin = polyy.min()
	inside_points = []
	for i in range(xmin,xmax+1):
		for j in range(ymin,ymax+1):
			if polypath.contains_point((i,j)) == True:
				inside_points.append((i,j))
	tock = time.time()
	#print "time for inside poly check is {0}".format(tock-tick)
	return inside_points

class hull_info:
	#points is a set of triples
	def __init__(self,points,param=300):
		self.points = points
		self.param = param
		###
		slice_dict = {}
		for triple in points:
			if triple[2] not in slice_dict.keys(): #intiate key
				#print 'key {} initated'.format(triple[2])
				slice_dict[triple[2]] = [(triple[0],triple[1])]
			
			elif triple[2] in slice_dict.keys(): #add to key if already exists
		
				old_values = slice_dict[triple[2]]
				old_values.append((triple[0],triple[1]))
				slice_dict[triple[2]] = old_values

		self.point_dict = slice_dict
		##
		filtered_dict={} #key is slices, returns filtered points. 
		vert_dict = {} #key is slices, returns vertices defing hull
		for k in slice_dict.keys():
			points = slice_dict[k]
			filtered = filter_out(points,param)
			filtered_dict[k] = filtered
			filtered = np.vstack(filtered) #change to matrix
			hull = ConvexHull(filtered)
			#plt.plot(filtered[:,0], filtered[:,1], 'o')
			#for simplex in hull.simplices:
			#	plt.plot(filtered[simplex, 0], filtered[simplex, 1], 'k-')
			vertices=[]
			for vert in hull.vertices:
				vertices.append(filtered[vert])
			vert_dict[k] = vertices #now returns list of vertices

		self.vert_dict = vert_dict
		self.filtered_dict = filtered_dict
		####

		inside_dict = {}
		for k in vert_dict.keys():
			inside = points_inside(vert_dict[k])
			inside_dict[k] = inside
		self.inside_dict = inside_dict #all points inside polygon fomres by vertices


	def density(self,slic):
		#will return the denisty of filtered points
		#inside polygon of slic
		return float(len(filtered_dict[slic]))/len(inside_dict[slic])
		
	def patch_together(self,STACK,save_dir=''):
		#will use hulls to patch together
		#using the CT stack, will use density
		#to decide what slice to use when hulls overlap
		#will save the image and return a dictionary telling you
		#what slice a pixel comes from
		Xs = np.shape(STACK)[0]
		Ys = np.shape(STACK)[1]
		patched = np.empty((Xs,Ys))
		patched[:] = -1
		best_density_dict={} #key will be coord, returns current density of hull taken from pixel
		what_slice = {} #will tell you what slice a point comes from
		for slic in self.inside_dict.keys(): #loop through slices
			for point in self.inside_dict[slic]: #through points in hull of that slice.
				if patched[point[0],point[1]] == -1:
					patched[point[0],point[1]] = STACK[point[0],point[1],slic]
					what_slice[point] = slic
					d = float(len(self.filtered_dict[slic]))/len(self.inside_dict[slic]) #density of hull
					best_density_dict[point] = d
				else:
					try:
						best = best_density_dict[point]
						new_d = float(len(self.filtered_dict[slic]))/len(self.inside_dict[slic])
						if new_d>best:
							patched[point[0],point[1]] = STACK[point[0],point[1],slic]
							best_density_dict[point] = new_d
							what_slice[point] = slic
						else:
							pass
					except KeyError:
						print 'key_error'
						pass
		#save
		image = Image.fromarray(patched)
		image.save('{0}match_with_patch_param={1}.tif'.format(save_dir,self.param),'tiff')
		return what_slice

		
		

	
	
def recover_homogenous_affine_transformation(p, p_prime):
	#Approximates affine transformation requited to get from
	#Vector p to p_prime

    # construct intermediate matrix
    Q       = p[1:]       - p[0]
    Q_prime = p_prime[1:] - p_prime[0]
    # calculate rotation matrix
    R = np.dot(np.linalg.inv(np.row_stack((Q, np.cross(*Q)))),
               np.row_stack((Q_prime, np.cross(*Q_prime))))

    # calculate translation vector
    t = p_prime[0] - np.dot(p[0], R)

    # calculate affine transformation matrix
    return R,t

def build_rot_M(theta):
	M = np.array(((np.cos(theta),-np.sin(theta)),(np.sin(theta),np.cos(theta))))
	return M


def approximate_affine(orig_list,list_onto):
	#takes two lists of tuples where we expect (list is the bsei points)
	#list is mapped on to list_onto with an affine transformation
	#and guesses the affine transformation by taking minimisation
	#over many points.
	def objective(AT):
		#AT = [l,theta,t1,t2]
		total = 0
		for i in range(len(orig_list)):
			x=np.array(orig_list[i])
			y=np.array(list_onto[i])
			total = total + np.linalg.norm(AT[0]*np.dot(build_rot_M(AT[1]),x)+np.array((AT[2],AT[3]))-y)

		return total

	x0 = [1,0,1,1]
	minimum = scipy.optimize.minimize(objective,x0,tol=0.01)	
	return minimum

def project_onto_ct(hull_info_ct,hull_info_bsei,chem_data,Aff_fiji):
		#takes hull_info classes of matched ct slice to the bsei image. Then projects the
		#chemical data onto the CT slice by using the Affine transformation from fiji (might have
		#to be inverse depending which way round it goes/
		#and BSEI hulls
		ct_in = hull_info_ct.inside_dict[1] #a list of tuples
		bsei_in = hull_info_bsei.inside_dict[1]
		chem_on_CT  = np.empty((stackX, stackY))## !!!!!stackX and stackY are globals define in main !!!!
		chem_on_CT[:] = -1

		opt = approximate_affine(hull_info_bsei.filtered_dict[1],hull_info_ct.filtered_dict[1])
		opt_param = opt.x
		print "opt_param is {0}".format(opt_param)
		for tup in bsei_in:
			tup = np.array(tup)
			#use my guess
			tup_prime = opt_param[0]*np.dot(build_rot_M(opt_param[1]),tup) + np.array((opt_param[2],opt_param[3]))

			#use fiji guess
			#tup_prime = np.dot(Aff_fiji[:,:2], tup) + Aff_fiji[:,2] #invert and +/- depends
			#on which way the fiji transformation was calcuated could be:
			#tup_prime = np.dot(Aff_fiji[:,:2], tup) + Aff_fiji[:,2]
		 
			chem_on_CT[int(round(tup_prime[0])),int(round(tup_prime[1]))] = chem_data[int(round
					(tup[0])),int(round(tup[1]))]
		
			
		ct_chem_image = Image.fromarray(chem_on_CT)
		ct_chem_image.save('CHEM_ON_CT.tif','tiff')
		ct_chem_image.close()
		
		
	
		


if __name__ == '__main__':
	##FOR PLANE FITTING##
	
	POINTS_DIR = './fiji_out_25_eye_four/'
	STACK_DIR =  '../by_eye_4_cropped/'
	SAVE_DIR = './result_by_eye_four/'
	STACK = DIL.load_image_sequence(STACK_DIR)
	ERROR = '25'
	SAVE_ALIGNED_STACK_DIR =SAVE_DIR+"/aligned_stack/"



	data = load_fiji_sift(POINTS_DIR)
	##put data as list of triples
	stackX,stackY,stackZ =np.shape(STACK)
	#plot_points(data,STACK,save_dir='./point_stack/')
	new_C = plane_fit(data,STACK,order=2,plot=False,save_name='{0}{1}_error'.format(SAVE_DIR,ERROR),save_stack=True,save_stack_dir=SAVE_ALIGNED_STACK_DIR) #works!
	#hull_class = hull_info(data)
	#hull_class.patch_together(STACK,SAVE_DIR)
	
	
	##FOR PROJECTING TO CT
	#CT_POINTS_DIR = './project_onto_ct/CT_points_for_subvolume1/'	
	#BSEI_POINTS_DIR = './project_onto_ct/BSEI_points_for_subvolume1/'
	
	#temp = Image.open('./CGA1-A_slice/match_to_scale_chem.tif')
	#chem = np.array(temp)
	#temp.close()
	
	
	#ct_points = load_fiji_sift(CT_POINTS_DIR)
	#bsei_points = load_fiji_sift(BSEI_POINTS_DIR)
	
	#ct_hull = hull_info(ct_points,param = 200000) # big param so doesnt filter at all
	#ct_hull.patch_together(STACK)
	#bsei_hull = hull_info(bsei_points,param =20000)

	

	#Aff_from_fiji = np.array(((0.99055301925573, 0.050353198184614, -113.8563174050355),(-0.070068216052093, 0.969428823434047, 158.65700844465664)))

	#project_onto_ct(ct_hull,bsei_hull,chem,Aff_from_fiji)


