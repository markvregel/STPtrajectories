def addhashtag(filename):
    inputfile = open(filename)
    outputfile = open(filename[:-2]+'_res.txt','w')
    for i in inputfile:
        outputfile.write("#'"+i)
    outputfile.close()
if __name__ == "__main__":        
    addhashtag('example_RTG.R')
    
