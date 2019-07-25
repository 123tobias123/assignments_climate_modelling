import statistics
import numpy

f=open("thedatafile.hhh", "r")
data = f.readlines()
intList=[]
floatList=[]
imaList=[]
imaListConverted = []
dictonaryList=[]
tuple_=()

# Read data into Lists/Dictonary
for i in data:
    line = i.rstrip()
    # Check at which section the loop is at
    if (i.startswith('#')):
        section = line # remembers the section that is being processed at the moment
        continue
    # Section 1: Integers
    if(section=='# section 1'):
        intList=line.split()
        
    # Section 2: Floats
    if(section=='# section 2'):
        floatList.append(line)

    # Section 3: Imaginary
    if(section=='# section 3'):
        imaList.append(line)
        if('j' not in line):
            tmp = line.split()
            imaListConverted.append(complex(float(tmp[0]),float(tmp[1])))
        else: 
            imaListConverted.append(complex(line))
    # Section 4: Dictonary
    if(section=='# section 4'):
        tmp = line.split()
        tuple_ = (tmp[0],tmp[1])
        dictonaryList.append(tuple_)

# Typecasting
intList = [int(x) for x in intList]
floatList = [float(x) for x in floatList]

# Output
# Section 1
print('Section 1')
print('mean: ', statistics.mean(intList))
print('variance: ',statistics.variance(intList))
print('standard deviation: '+str(statistics.stdev(intList)))
# Section 2
print('Section 2')
print('mean: '+str(statistics.mean(floatList)))
print('variance: '+str(statistics.variance(floatList)))
print('standard deviation: '+str(statistics.stdev(floatList)))
# Section 3
print('Section 3')
print('Product of complex numbers: ' + str(numpy.prod(imaListConverted))) 
# Section 4
print('Section 4')
[print(str(x[0])+ ' ' + str(x[1])) for x in sorted(dictonaryList)]

f.close()
"""
TFL: I really like this. It's very well organised and clean code. I only have a few very minor comments. 
 1) Use white space consistently (lines 4-11). But this is a personal preference. :)
 2) You used multiple if statements in your reading loop to check the current section. This is needed for the first and second if blocks, but it's slightly inefficient for the others because it will check every section on every iteration of the loop. If one uses elif, then once a tested expression is true, the following elif/else expressions are ignored. 
 3) You didn't use a dictionary data type for section 4; it would have made sense to do so. However, it's perfectly fine as is, as a dictionary ultimately isn't necessary for the task.
 4) You can also pass an arbitrary number of arguments to print(), and they need not be strings. eg. print('mean: ', statistics.mean(floatList))
 5) It could be good idea to close f when it's not needed, so that other processes can use it.
"""