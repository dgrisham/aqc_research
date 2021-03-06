#!/usr/bin/env python2
# author: mberntso  last update: 3/28/16

# current input is for sample data as supplied by Prof. Coffey  

# this algorithm could be improved by implementing a file format
# for data input

from operator import sub

# initialize all tasks as scheduled

## numTasks = 6
x = [1,1,1,1,1,1]

## global indeces for first and last scheduled tasks
firstIndex = 0
lastIndex = len(x)-1

# initialize start times and precedence matrix later

# model inputs

# duration

d = [5,5,5,5,5,5]

## slew time in accordance with index ordering
s = [[0 for z in range(len(x))] for z in range(len(x))]

### s[i][j] = time from i to j
s[0][1] = 3
s[0][2] = 4
s[0][3] = 1
s[0][4] = 7
s[0][5] = 10
s[1][2] = 5
s[1][3] = 4
s[1][4] = 2
s[1][5] = 5
s[2][3] = 5
s[2][4] = 4
s[2][5] = 8
s[3][5] = 4
s[3][5] = 1
s[4][5] = 3
# for i in range(len(s)-1):
#     for j in range(i+1, len(s)):
#         s[i][j] = 1

# window begin and end times
wB = [103,105,300,100,200,300]
wE = [110,111,600,1200,400,350]

# priorities
p = [1,2,3,1,2,3]

# very large value (not used currently)
M = 999

# build y (precedence) matrix (not used currently, b/c ordering is assumed)
def precedenceMatrix(x):
    numTasks = len(x)

    y = [[0 for z in range(numTasks)] for z in range(numTasks)]

    # y[i][j] = 1 if i < j

    # todo: for loops, when data from user input

    for i in range(len(y)-1):
        for j in range(i+1, len(y)):
            y[i][j] = 1

    return y

# build bounds for starting time b
bMin = wB
bMax = map(sub, wE, d)

# findNext/Prev scheduled task (index)
def findNext(x, index, last):
    if index >= last:
       return last
    elif x[index+1] == 1:
        return index+1
    else:
        return findNext(x, index+1, last)

def findPrev(x, index, first):
    if index <= first:
        return first
    elif x[index-1] == 1:
        return index-1
    else:
        return findPrev(x, index-1, first)

# remove task at (index) and update bounds on indeces
def removeTask(x, index, first, last, currentV):
    x[index] = 0
    currentV = currentV - p[index]
    if index == first:
        first = findNext(x, index, last)
    if index == last:
        last = findPrev(x, index, first)
    return (x,first,last,currentV)

# intialize the objective function V as if all tasks are scheduled
currentV = sum(p)


# recursively determine the optimal schedule
#
# Perhaps following an approach in line with the additive algorithm
# would clean up the mess of conditionals in this function
# Update (4/12/16): Additive Algorithm is just as complex and problem dependent
#
# See Figure 1 and 1a in the Geoffrion paper for a similar algorithm.
# The Geoffrion paper also includes some explanation of a nonlinear
# objective function for binary free variables.  That seems like a plausible
# extension to what is here.
def optimalValue(x, index, first, last, currentV, bMin, bMax):
    while index <= last: # O(n)
        next = findNext(x, index, last) # O(n)
        if index == last:
            if bMin[index] > bMax[index]:
                (x,first,last,currentV) = removeTask(x, index, first, last, currentV)
            break
        if bMin[index] > bMax[index]:
            (x,first,last,currentV) = removeTask(x, index, first, last, currentV)
        else:
            if bMin[index]+d[index]+s[index][next] > bMax[next]:

                # removal and comparison

                # create copies of bMin and bMax
                bMin_Z = list(bMin)
                bMax_Z = list(bMax)

                (x_A,first_A,last_A,currentV_A) = removeTask(list(x), index, first, last, currentV)
                index_A = findPrev(x, index, first)
                # store result of removing current task as B
                (x_B,first_B,last_B,currentV_B,bMin_B,bMax_B) =\
                optimalValue(x_A, index_A, first_A, last_A, currentV_A, bMin_Z, bMax_Z) # recursive => O(n^2)

                (x_C,first_C,last_C,currentV_C) = removeTask(list(x), next, first, last, currentV)
                index_C = findNext(x, index, last)
                # store result of removing next task as D
                (x_D,first_D,last_D,currentV_D,bMin_D,bMax_D) =\
                optimalValue(x_C, index_C, first_C, last_C, currentV_C, bMin_Z, bMax_Z)

                if currentV_B >= currentV_D:
                    if currentV_B == currentV_D:
                        if p[index] >= p[next]:        
                            (x,first,last,currentV, bMin, bMax) = (x_B,first_B,last_B,currentV_B,bMin_B,bMax_B)
                            index = findPrev(x, index, first)
                        else:
                            (x,first,last,currentV, bMin, bMax) = (x_D,first_D,last_D,currentV_D,bMin_D,bMax_D)
                    else:
                        (x,first,last,currentV, bMin, bMax) = (x_B,first_B,last_B,currentV_B,bMin_B,bMax_B)
                        index = findPrev(x, index, first)
                else:    
                    (x,first,last,currentV, bMin, bMax) = (x_D,first_D,last_D,currentV_D,bMin_D,bMax_D)

            else:

                # adjust bounds if infeasible overlap
                if bMin[index]+d[index]+s[index][next] <= bMin[next]:
                    if bMax[index]+d[index]+s[index][next] > bMin[next]:
                        if bMax[index]+d[index]+s[index][next] > bMax[next]:
                            bMax[index] = bMax[next]-(d[index]+s[index][next])
                else:
                    bMin[next] = bMin[index]+(d[index]+s[index][next])

        index = findNext(x, index, last)

    return (x,first,last,currentV, bMin, bMax)


# print the results
(sched, blah, bleh, objFunc, startMin, startMax) = optimalValue(x, firstIndex, firstIndex, lastIndex, currentV, bMin, bMax)

print('Scheduled Tasks: ')
print(sched)
print('Objective Function Value: ')
print(sum([a*b for a,b in zip(sched,p)]))
print('Start time min bound: ')
print(startMin)
print('Start time max bound: ')
print(startMax)






