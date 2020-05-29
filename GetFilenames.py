

#get filenames
fnames=[]
import os
TargetPath=""#"/Users/cameronerickson/Desktop/Academics/NPL/dePol/Data Analysis/NMR_Data/12-18-2017-Mon/"
Target="../T1Measurement3"#75Burst2Vpp"#Background75Burst#
for file in os.listdir(TargetPath+Target):
    if file.endswith(".txt"):
        print(os.path.join(file))
        fnames.append(os.path.join(file))

#print(int(fnames[0].replace("-","")[:-4]))

#Sort the filenames from oldest to newest
fnames =sorted(fnames, key=lambda fname: int(fname.replace("-","")[:-7]))
#print(fnames)

#Creating File
path=""#"""/Users/cameronerickson/Desktop/Academics/NPL/dePol/Data Analysis/NMR_Data/12-15-2017-Fri/"
fname="FileNames.txt"


f=open(path+fname,'w')

for line in fnames:
    f.write("../"+TargetPath+Target+"/"+line)
    f.write('\n')

f.close()
