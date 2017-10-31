# import os

# directory_to_check = "/media/ravindra/SAMSUNG/Rinkle/Ge-based/As2GeCd" # Which directory do you want to start with?

# def my_function(directory):
#       #print("Listing: " + directory[31:])
#       print("\t-" + "\n\t-".join(os.listdir("."))) # List current working directory
#       #os.system('python bandplotting-orb-resolved-bug-removed.py ')

# # Get all the subdirectories of directory_to_check recursively and store them in a list:
# directories = [os.path.abspath(x[0]) for x in os.walk(directory_to_check)]
# #directories.remove(os.path.abspath(directory_to_check)) # If you don't want your main directory included

# for i in directories:
#       os.chdir(i)         # Change working Directory
#       my_function(i)      # Run your function



import glob, os

def scanfolder():
    for path, dirs, files in os.walk('/home/ravindra/Rinkle-pressure-bands-with-soc-primitive'):
        #print path
        #print dirs
        #print files
        for _dir in dirs:
        	os.chdir(path+"/"+_dir)
        	if os.path.isfile('vasprun.xml') :
        		os.system('cp /home/ravindra/Rinkle-pressure-bands-with-soc-primitive/bandplotting-orb-resolved-bug-removed-without-borders.py .' )
        		cwdir1 = os.getcwd().replace("/","/")
        		#print path
        		print "plotting " + path + cwdir1[99:]
        		cwdir2 = os.getcwd().replace("/","-")
        		#print cwdir2[90:]
        		os.system("time python bandplotting-orb-resolved-bug-removed-without-borders.py " + cwdir2[90:])

scanfolder()