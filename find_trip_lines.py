import csv
import numpy as np
name = '3D_block.csv'
#np.set_printoptions(threshold=np.inf)
#Given a file name, return an array of parameter "value"
def get_value_array(name,value):

	with open(name, 'r') as csvfile:
		reader = csv.reader(csvfile)	
		rownum = 0
		max_x,max_y,max_z = 0,0,0
		data_array = []
		for row in reader: #iterate through each row
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

#Find a set that represents the 3 neighbors of a point on a TJ
#The set contains the FID of the nearest neighbors 

def find_neighbor_set(i,j,k,FID_array):
	size_z,size_y,size_x = FID_array.shape
	neighbor_set = set()
	for zpos in xrange(i - 1,i + 2,1):
		for ypos in xrange(j - 1,j + 2, 1):
			for xpos in xrange(k - 1,k + 2, 1):
				neighbor = FID_array[zpos % size_z][ypos % size_y][xpos % size_x]
				neighbor_set.add(neighbor)
	return neighbor_set
				
#returns (index of list_of_NN in neighbor list,1 = found in NN list 0 = not found, 
#append onto end of NN list)

def find_index(neighbor_set, list_of_NN):
	for i in xrange(len(list_of_NN)):
		if (list_of_NN[i] == neighbor_set):
			return (i,1)
	return (len(list_of_NN),0) #Couldn't find in NN list, so append onto end

#label_TJ takes in the triple junction array 
#and the feature ID array. It uniquely identifies each TJ line 
#by it's 3 neighbors. 

#TJ_Num_array --> -1 = Not a TJ point
# Integer greater than or equal to zero = TJ number of point


def label_TJ(TJ_array, FID_array):
	size_z,size_y,size_x = TJ_array.shape
	TJ_Num_array = np.ones(((size_z,size_y,size_x)))*-1 #-1 denotes not a triple junction point
	TJ_Num_array_linear = np.zeros((size_z*size_y*size_x,4)) #TJ_Num_array in form [ [x,y,z, TJ_Num], ...] 
	list_of_NN = np.array([])# List of sets of triple junction tri-neighbor combinations
	count = 0
	for i in xrange(size_z):
		for j in xrange(size_y):
			for k in xrange(size_x):
				TJ_Num_array_linear[count] = [k,j,i,-1] #Add current point to linear array.
				if (TJ_array[i][j][k] == 0): #Found a triple junction point!
					neighbor_set = find_neighbor_set(i,j,k,FID_array) #get set of neighbors to this point
					list_index = find_index(neighbor_set,list_of_NN) #get index of these neighbors in list_of_NN		
					TJ_Num_array[i][j][k] = list_index[0] #Set this point equal to that index
					if (list_index[1] == 0): #This neighbor_set was not in list_of_NN
						list_of_NN = np.append(list_of_NN,neighbor_set) #So we have to append it to the list of NN
					TJ_Num_array_linear[count][3] = list_index[0]#update linear array
				count += 1
	return TJ_Num_array_linear

value = "TJEuclideanDistances"
TJ_array = get_value_array(name,value)
print TJ_array.shape
value = "FeatureIds"
FID_array = get_value_array(name,value)
TJ_Num_array_linear = label_TJ(TJ_array, FID_array)
np.savetxt('3D_Block_TJ_lines.out',TJ_Num_array_linear,delimiter=',', fmt = '%d')

