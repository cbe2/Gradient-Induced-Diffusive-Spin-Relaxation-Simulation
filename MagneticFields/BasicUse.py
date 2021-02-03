import magpylib as magpy

# Millimeter for lengths
# Degree for angles
# Millitesla for magnetization/remanence, magnetic moment and magnetic field,
# Ampere for currents.

#create sources
#magnetization is in x dir with magnitude 1 in CHECK UNITS LATER
#dimesions are [diameter, height]=[4,5] (mm)
#make cylinder with center 5 mm above origin in z direction
s1 = magpy.source.magnet.Cylinder( mag = [1,0,0],dim = [10,5], pos = [0,0,5])

#makes loop in x-y plane
s3 = magpy.source.current.Circular( curr = 1, dim =10) #current is one amp clockwise (in z direction), dimension is 10 mm

#create collection
c = magpy.Collection(s1,s3) #now the system can be moved and manipulated as one

#display system

#specify markers to show places of interest:
markerPos = [(0,0,0,'origin'),(10,10,10),(-10,-10,-10)]

#Evaluate field at orign
print(c.getB([0,0,0]))

#show the collection
magpy.displaySystem(c,markers=markerPos,direc=True)
