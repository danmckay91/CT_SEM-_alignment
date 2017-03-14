from PIL import Image
import os, glob
import numpy as np
import os.path
import sys
import math
import itertools
import time



def load_image_sequence(dir):
	stack_string=[]
	for tif in glob.glob("{0}*.tif".format(dir)):
		stack_string.append(tif)
	
	stack_string.sort() 
	
	stack_list=[]
	for tif in stack_string: #loads each tif in stack as a numpy array of grey values
		temp = Image.open(tif)
		stack_list.append(np.array(temp))
		temp.close()
	stack=np.dstack(stack_list)
	return stack
	

def save_as_im_sequence(a_stack,dir='./'):
	for i in range(np.shape(a_stack)[2]):
		slice = a_stack[:,:,i]
		result = Image.fromarray(slice)
		result.save('{0}/{1}.tif'.format(dir,i))
		
		
def return_slice_no(astack,theta,phi,kstar):
	#reconstructs a slice kstar after rotation theta about x-axis
	#and phi about y-axis
	#star denotes a rotated coordiante
	xdim, ydim, zdim =np.shape(astack)
	icent,jcent,kcent = ((xdim+1)/2, (ydim+1)/2, (zdim+1)/2)
	slice = np.zeros((xdim,ydim))
	RotM =np.array([[np.cos(phi), 0, np.sin(phi)],[np.sin(theta)*np.sin(phi),np.cos(theta),-np.sin(theta)*np.cos(phi)],[-np.cos(theta)*np.sin(phi), np.sin(theta), np.cos(theta)*np.cos(phi)]]) 
	RotMinv = np.linalg.inv(RotM)
	for istar in range(xdim):
		for jstar in range(ydim):
			[i,j,k] = np.dot(RotMinv,np.array([istar-icent,jstar-jcent,kstar-kcent]))+np.array([icent,jcent,kcent])
			i=int(round(i))
			j=int(round(j))
			k=int(round(k))
			if i<ydim and i>=0 and j<xdim and j>=0 and k<zdim and k>=0: #make sure we are not out of bounds 
				slice[istar,jstar] = astack[i,j,k] #round is nearest neighbour
			else:
				pass
	return slice
	
def rotate_stack(astack,RotM):
	RotMinv = np.linalg.inv(RotM)
	xdim, ydim, zdim =np.shape(astack)
	icent,jcent,kcent = ((xdim+1)/2, (ydim+1)/2, (zdim+1)/2)
	max_dim = max([xdim,ydim,zdim])
	print max_dim
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
	
			
	
	# for i in range(xdim):
		# for j in range(ydim): #the + centre[2] might not be needed due to slice numbering. i.e rotating about the origin!!
			# newk = k/(math.cos(theta)*math.cos(phi)) + (i-centre[0])*math.tan(theta)/(math.cos(phi))  - math.tan(theta)*(j-centre[1])/math.cos(phi)# retruns slice number for each i,j after rotation

			# if r <zdim and r>=0:
				# slice[i,j] = astack[i,j,int(r)] # fills in new slice
			# else: 
				# pass
	# return slice	

def crop(image,x,y,i,j):
	'''crops an image to rectangle x,y with top left corner
	corner (i,j)	
	'''
	#output=np.zeros((x,y))
	output = image[i:x+i,j:y+j]
	return output	
	
def global_pc(im1,im2):
	a= im1.tolist()
	b = im2.tolist()
	return  np.cov(a,b)[0,1]
	
	
def find_coords(a_stack):
	sz=np.shape(a_stack)
	max_area=int(math.ceil(max(math.sqrt(sz[0]**2+sz[0]**2),math.sqrt(sz[1]**2+sz[1]**2),math.sqrt(sz[2]**2+sz[2]**2),math.sqrt(sz[0]**2+sz[1]**2),math.sqrt(sz[0]**2+sz[2]**2),math.sqrt(sz[1]**2+sz[2]**2))))
	nd1,nd2,nd3 = max_area,max_area,max_area
	product=[]
	for item in itertools.product(range(nd1),range(nd2),range(nd3)):
		product.append(item)
		
	product=np.array(product)
	return product	

	
def reslice_stack(a_stack,rotM,COR,o=0,pad=True):
	'''will reslice a_stack after rotating theta about axis ax1, phi about ax2 and alpha about ax3.
	method is 0, 1 ,2 ,3 splice interpolation
	'''
	a_stack=np.uint8(a_stack)
	product= find_coords(a_stack)
	tick=time.time()

	if pad==True: #imbeds smaller stack in larger one 
		s = max(sz)
		max_area=int(math.ceil(max(math.sqrt(sz[0]**2+sz[0]**2),math.sqrt(sz[1]**2+sz[1]**2),math.sqrt(sz[2]**2+sz[2]**2),math.sqrt(sz[0]**2+sz[1]**2),math.sqrt(sz[0]**2+sz[2]**2),math.sqrt(sz[1]**2+sz[2]**2))))
		imagepad = np.zeros((max_area,max_area,max_area))
		from_edge=[math.floor((max_area - l)/2) for l in sz]
		imagepad[from_edge[0]:from_edge[0]+sz[0],
		from_edge[1]:from_edge[1]+sz[1],
		from_edge[2]:from_edge[2]+sz[2]] = a_stack 
	else:
		imagepad=a_stack
		
	
	
	nd1,nd2,nd3 = np.shape(imagepad) #find middles
	midx = COR[0]
	midy = COR[1]
	midz = COR[2]
	mid=np.array([midx,midy,midz])
	##
	'''
	tick = time.time()
	product=[]
	for item in itertools.product(range(nd1),range(nd2),range(nd3)):
		product.append(item)
		
	product=np.array(product)
	'''
	###
	rotated_coord=np.dot(product-mid,rotM)+mid
	
	
	xout = rotated_coord[:,0]
	yout = rotated_coord[:,1]
	zout = rotated_coord[:,2]
	

	rotated=map_coordinates(imagepad,[xout,yout,zout],order=o)
	tock=time.time()
	print 'Rotation complete. total time took {0}'.format(tock-tick) 
	return np.uint8(np.reshape(rotated,np.shape(imagepad)))
		
	
if __name__ == "__main__":	
	stack=load_image_sequence('./test_stack/')
	# rotated = np.zeros(np.shape(stack))
	# for kstar in range(np.shape(stack)[2]):
		# rotated[:,:,kstar] = return_slice_no(stack,np.pi/2,0,kstar)
	# save_as_im_sequence(rotated,'./rotated/')
	slice = return_slice_no(stack,np.pi/180*5,np.pi/180*5,100)
	result =Image.fromarray(slice)
	result.save('./test_slice/match_to.tif', 'tiff')



		
#(-0.004362138360144,0.0200999785901,47,22,146)		
		
		
		
		
