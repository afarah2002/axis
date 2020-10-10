import glob
import os
import shutil

version = 1
src = "/media/nasa01/5A8C680F8C67E3CB/Documents and Settings/Albert/My Documents/Howe Internship/firestorm/main/current ion thruster" + "/V" + str(version) + "_stls/"

# make new director next to converter file
stl_directory = "V" + str(version) + "_gdmls/"
parent_dir = "/home/nasa01/Documents/CAD2GEANT4Converter/"
stl_dest = parent_dir + stl_directory
os.system("mkdir " + stl_dest)
print("Directory '% s' created" % stl_dest) 

# copy files to converter directory
os.system("cp '" + src + "'*.STL " + stl_dest)

# # make GDML directory
# gdml_directory = "V" + str(version) + "_gdmls/"
# os.system("mkdir " + parent_dir + gdml_directory)

# assign materials to stls
stl_files = glob.glob(stl_dest + "*.STL")

material = "_Aluminum"
for src in os.listdir(stl_dest):
	src1 = str(src).replace(" ", "")
	src1 = str(src1).replace("_", "")
	src1 = str(src1).replace("-", "")
	dst = src1[0:-4] + material + src1[-4:]
	os.rename(stl_dest + src, stl_dest + dst)

# convert to gdml
converter_path = parent_dir + "stl_gdml.py "
axis_dest = "/home/nasa01/Documents/howe_internship/axis/"
for src in os.listdir(stl_dest):
	print src
	# dst = axis_dest 
	name = src[0:-4] + " "
	# print(dst)
	os.system("python3 " + converter_path + name + stl_dest + src)
	os.system("cp " + stl_dest + "* " + axis_dest)

# copy gdmls to axis directory
# os.system("mkdir " + axis_dest)
# os.system("cp " + stl_dest + "*.gdml " + axis_dest)