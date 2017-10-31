import glob, os

def scanfolder():
    for path, dirs, files in os.walk('/media/ravindra/SAMSUNG/Rinkle'):
        #print path
        #print dirs
        #print files
        for _dir in dirs:
        	os.chdir(path+"/"+_dir)
        	if os.path.isfile('vasprun.xml'):
        		#os.system('cp /media/ravindra/SAMSUNG/Rinkle/bandplotting-orb-resolved-bug-removed-without-borders.py .' )#+ " -a "+ faa_file_path + " -b " + fna_file_path)
        		cwdir = os.getcwd().replace("/","-")
        		#print path + cwdir[31:]
        		pdffilenamebig = str(glob.glob('*.pdf'))
        		print cwdir[31:]+'.pdf'
        		os.system("cp " + cwdir[31:]+'-reduced.pdf  '+ '/home/ravindra/pending_papers/manuscript-intermetallics/SI/')
        		#os.system("ps2pdf " + cwdir[31:]+'.ps  ' + cwdir[31:]+'-reduced.pdf')
        		#os.system("pdf2ps " + cwdir[31:]+'.pdf')

scanfolder()