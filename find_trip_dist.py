import csv
import numpy as np
from evtk.hl import pointsToVTK
from evtk.hl import imageToVTK 

name = '3D_Block_complex.csv'
output = './3D_Block_complex_2'

#np.set_printoptions(threshold=np.inf)
#Given a file name, return an array of parameter "value"
#in array[z][y][x] form
#Assume output csv format from paraviewi
def get_value_array(name,value):

	with open(name, 'r') as csvfile:
		reader = csv.reader(csvfile)	
		rownum = 0
		max_x,max_y,max_z = 0,0,0
		data_array = []
		for row in reader:
			if (not rownum): #Identify columns via first row
				V_col = row.index(value)
				X_col = row.index('Structured Coordinates:0')
				Y_col = row.index('Structured Coordinates:1')
				Z_col = row.index('Structured Coordinates:2')
			else:
				xval,yval,zval = int(row[X_col]),int(row[Y_col]),int(row[Z_col])
				if xval > max_x: max_x = xval
				if yval > max_y: max_y = yval
				if zval > max_z: max_z = zval
				data_array += [row]
			rownum += 1
		
		V_array = np.ones(((max_z + 1,max_y + 1,max_x + 1)))*-1
		for row in data_array:
			xval,yval,zval,V = int(row[X_col]),int(row[Y_col]),int(row[Z_col]),float(row[V_col])
			V_array[zval][yval][xval] = V
			
	return V_array

#Given a point p (z,y,x), find it's 
#root in the union array

def find_root(p, union_array):
	while (make_char_num(p) != make_char_num(union_array[p[0]][p[1]][p[2]]) ):
		p = (union_array[p[0]][p[1]][p[2]])
	return p

#Given two points, p1 and p2, in the form (z,y,x)
#and the union_array encoding the root of every
#point, union the two points.

def union(p1,p2,union_array):
	root1 = find_root(p1,union_array)
	root2 = find_root(p2,union_array)
	p1_v,p2_v = make_char_num(root1),make_char_num(root2)
	if (p1_v < p2_v):
		union_array[root2[0]][root2[1]][root2[2]] = root1
	elif (p2_v < p1_v):
		union_array[root1[0]][root1[1]][root1[2]] = root2
	return union_array

#Check whether val1 and val2, in (z,y,x) form, 
#have the same feature ID
def same_FID(val1,val2,FID_array):
	FID1 = FID_array[val1[0]][val1[1]][val1[2]]
	FID2 = FID_array[val2[0]][val2[1]][val2[2]]
	return (FID1 == FID2)

#Makes a number to represent a point (z,y,x)
#The "characteristic number" of the point
#Assumes x,y,z < 1000
def make_char_num(p):
	(z,y,x) = p
	return z*(10**6) + y*(10**3) + x

#If (z,y,x) is on grain boundary/ triple point, ect.
#coord = (z,y,x). 0 = Not a special point. 1 = is a special point.
def is_special_point(array,coord):
	val = array[coord[0]][coord[1]][coord[2]][0]
	return (val + 1) 

"""Build a 3D union data structure of material 
Tuple at every point. (-1,-1,-1) if not on a GB
(x,y,z) of root of point if on a grain boundary"""


def Union_Array(TJ_array,FID_array):
	size_z,size_y,size_x = TJ_array.shape
	union_array = np.ones( ((size_z,size_y,size_x)), dtype = ('i4,i4,i4') )
	union_array.fill( (-1,-1,-1) ) 
	#label each TJ point on the material with a unique number
	for z in range(size_z):
		for y in range(size_y):
			for x in range(size_x): #Label all points on GB as non-zero and distinct in union_array
				if (not TJ_array[z][y][x]): #Distance is 0 from a GB, so we're on one
					union_array[z][y][x] = (z,y,x) #Root is current position 

	#(-1,-1,-1) denotes a non-GB point. Else, (z,y,x) denotes root of current point	
	for z in range(size_z):
		for y in range(size_y):
			for x in range(size_x):
				val = union_array[z][y][x]
				if (is_special_point(union_array,(z,y,x))): #On a triple junction point	
					union_array = union_point_with_neighbors((z,y,x),union_array,FID_array)
	return union_array

#Given a point p(z,y,x) and the union_array,
#return a union_array with p union'd with 
#'special point' neighbors

def union_point_with_neighbors(p,union_array,FID_array):
	z,y,x = p
	size_z,size_y,size_x = union_array.shape
	val = union_array[z][y][x]
	for z_step in range(-1,2,1):
		for y_step in range(-1,2,1):
			for x_step in range(-1,2,1):
				pid = ( (z + z_step) % size_z ,(y + y_step) % size_y , (x + x_step) % size_x )
				if ( (is_special_point(union_array,pid)) and (same_FID(val,pid,FID_array)) ):
					union_array = union(val,pid,union_array)
	return union_array

#Given a union array, make a 1d list of (x,y,z) 
#points and their grain boundary/ triple point
#number

def make_TJ_list(union_array,TJ_array):
	size_z,size_y,size_x = TJ_array.shape
	TJ_list = []
	for z in range(size_z):
		for y in range(size_y):
			for x in range(size_x):
				p = (z,y,x)
				if (is_special_point(union_array, p)):
					root = find_root(p,union_array)
					root_num = make_char_num(root)
					TJ_list += [[z,y,x,root_num]]	
				else:
					TJ_list += [[z,y,x,-1]]
	TJ_list = sorted(TJ_list,key=lambda x: x[3])
	return TJ_list


#Take an ordered list, and make all numbers consecutive integers
def clean_list(L):
	count = -1
	last_num = L[0][3] 
	L[0][3] = count
	for i in range(1,len(L)):  
		curr_num = L[i][3] #current value 

		if (curr_num > last_num):
			count += 1
			L[i][3] = count

		else: L[i][3] = count

		last_num = curr_num

	return np.array(L)
 

def output_to_vtk_points(L):
	z = np.array(L[:,0])
	y = np.array(L[:,1])
	x = np.array(L[:,2])
	vals = np.array(L[:,3])
	
	pointsToVTK(output,x,y,z, data = {"TJ_dist" : vals})

def output_to_vtk_image(L,size_x,size_y,size_z):
	vals = np.zeros(((size_x,size_y,size_z)))
	for i in xrange(len(L)):
		x,y,z = L[i][2],L[i][1],L[i][0]
		val = L[i][3]
		vals[x][y][z] = val
	imageToVTK(output, cellData = {"TJ Num" : vals})	
		

value = "TJEuclideanDistances"
TJ_array = get_value_array(name,value)
value = "FeatureIds"
FID_array = get_value_array(name,value)

size_z,size_y,size_x = TJ_array.shape
union_array = Union_Array(TJ_array,FID_array)
TJ_list = make_TJ_list(union_array,TJ_array)
TJ_list = clean_list(TJ_list)
print TJ_list
output_to_vtk_image(TJ_list,size_x,size_y,size_z)
#output_to_vtk(TJ_list)
