import csv
import numpy as np
from evtk.hl import pointsToVTK
name = 'test2.csv'
output = '2d_GB_list_test.csv'

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


def union_array(GB_array,FID_array):
	size_z,size_y,size_x = GB_array.shape
	union_array = np.ones( ((size_z,size_y,size_x)), dtype = ('i4,i4,i4') )
	union_array.fill( (-1,-1,-1) ) 
	#label each GB point on the material with a unique number
	for z in range(size_z):
		for y in range(size_y):
			for x in range(size_x): #Label all points on GB as non-zero and distinct in union_array
				if (not GB_array[z][y][x]): #Distance is 0 from a GB, so we're on one
					union_array[z][y][x] = (z,y,x) #Root is current position 

	#(-1,-1,-1) denotes a non-GB point. Else, (z,y,x) denotes root of current point	
	for z in range(size_z):
		for y in range(size_y):
			for x in range(size_x):
				val = union_array[z][y][x]
				if (is_special_point(union_array,(z,y,x))): #On a triple junction point	
					rid = (z,y,(x + 1) % size_x) #index of point to right
					did = (z, (y + 1) % size_y,x) #index of down point 
					drid = (z, (y + 1) % size_y,(x + 1) % size_x) #index of down-right point
					urid = (z, (y - 1) % size_y,(x + 1) % size_x) #index of up-right point 
					if (is_special_point(union_array,rid)) and (same_FID(val,rid,FID_array)): union_array = union(val,rid,union_array)
					if (is_special_point(union_array,did)) and (same_FID(val,did,FID_array)): union_array = union(val,did, union_array)
					if (is_special_point(union_array,drid)) and (same_FID(val,drid,FID_array)): union_array = union(val, drid, union_array)
					if (is_special_point(union_array,urid)) and (same_FID(val,urid,FID_array)): union_array = union(val, urid, union_array)
	return union_array


#Given a union array, make a 1d list of (x,y,z) 
#points and their grain boundary/ triple point
#number

def make_GB_list(union_array):
	size_z,size_y,size_x = GB_array.shape
	GB_list = []
	for z in range(size_z):
		for y in range(size_y):
			for x in range(size_x):
				p = (z,y,x)
				if (is_special_point(union_array, p)):
					root = find_root(p,union_array)
					root_num = make_char_num(root)
					GB_list += [[z,y,x,root_num]]	
				else:
					GB_list += [[z,y,x,-1]]
	GB_list = sorted(GB_list,key=lambda x: x[3])
	return GB_list


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
 

def output_to_vtk(GB_list):
	z = np.array(GB_list[:,0])
	y = np.array(GB_list[:,1])
	x = np.array(GB_list[:,2])
	vals = np.array(GB_list[:,3])
	
	pointsToVTK("./points",x,y,z, data = {"GB_dist" : vals})

value = "GBEuclideanDistances"
GB_array = get_value_array(name,value)
value = "FeatureIds"
FID_array = get_value_array(name,value)

size_z,size_y,size_x = GB_array.shape
union_array = union_array(GB_array,FID_array)
GB_list = make_GB_list(union_array)
GB_list = clean_list(GB_list)

"""
f = open(output,'wt')
writer = csv.writer(f)
writer.writerow( ('x','y','z','num'))
for i in range(len(GB_list)):
	writer.writerow( (GB_list[i][0],GB_list[i][1],GB_list[i][2],GB_list[i][3] ))

"""
output_to_vtk(GB_list)
np.savetxt('2d_GB_list_test.txt', GB_list, delimiter=',', fmt = '%.4f')
