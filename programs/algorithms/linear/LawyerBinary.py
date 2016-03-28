#!/usr/bin/env python2
# author: mberntso  last update: 1/23/15  

class Cell:
    def __init__(self, assigned, lawyer, case, time):
        self.assigned = assigned
        self.lawyer = lawyer
        self.case = case
        self.time = time
        return


def organizeData(data):

    table = [[Cell(0,i,j,data[i][j]) for j in range(len(data[0]))] for i in range(len(data))]
    #for cellset in table:
    #    for cell in cellset:
    #        print("lawyer: {}".format(cell.lawyer))
    #        print("case: {}".format(cell.case))
    #        print("time: {}".format(cell.time))

    # data table and lists by lawyer created

    # now create sorted list of times in ascending order for additive algorithm

    sortedList = []

    for i in range(len(data)):
        for j in range(len(data[0])):
            sortedList.append(table[i][j])
    
    sortedList = sorted(sortedList, key=lambda Cell: Cell.time)
            
    #for i in range(len(sortedList)):
    #    print(sortedList[i].time)

    return table, sortedList


def updateIncumbent(table, alloc, incumbent):

    # check if updated allocation that satisfies constraints gives better obj func value

    newValue = 0

    for i in range(len(alloc)):
        newValue += table[i][alloc[i]].time

    print(newValue)

    if newValue < incumbent:
        incumbent = newValue

    return incumbent


def updateAlloc(table):

    # update alloc from cells in table

    alloc = [-1]*len(timeData)
    for i in range(len(table)):
        for j in range(len(table[0])):
            if table[i][j].assigned == 1:
                alloc[i] = j

    print(str(alloc) + 'alloc')

    return alloc

# recursively determine the optimal allocation of cases
# using Balas additive algorithm, assume root node cannot be solution

def optimalAlloc(table, sortList, alloc, incumbent, index):
    
    # check constraints (one case per lawyer)

    sumList = []

    for i in range(len(table)):
        tempSum = 0
        for j in range(len(table[0])):
            tempSum += table[i][j].assigned

        if tempSum > 1:  # lawyer assigned to an impossible number of cases
            print('impossible')
            return (alloc, incumbent)

        sumList.append(tempSum)

    print(str(sumList)+'list')
    for i in sumList:
        # Check for infeasible constraint

        if i == 0:
            # Depth First Recursion
            table[sortList[index].lawyer][sortList[index].case].assigned = 1
            alloc = updateAlloc(table)
            print('A')
            (allocA, incumbentA) = optimalAlloc(table, sortList, alloc, incumbent, index+1)
            print(allocA, incumbentA)

            table[sortList[index].lawyer][sortList[index].case].assigned = 0
            table[sortList[index+1].lawyer][sortList[index+1].case].assigned = 1
            alloc = updateAlloc(table)
            print('B')
            (allocB, incumbentB) = optimalAlloc(table, sortList, alloc, incumbent, index+1)
            print(allocB, incumbentB)

            if incumbentB < incumbentA:
                return (allocB, incumbentB)
            return (allocA, incumbentA)

    # satisfied constraints
    alloc = updateAlloc(table)
    incumbent = updateIncumbent(table, alloc, incumbent)

    return (alloc, incumbent)


# Results:

# Setup data:

# table of time data for each case (col) by each lawyer (row) => data[row][col]

timeData = [[1., 2., 3.],
        [4., 5., 6.],
        [7., 8., 9.]]




#[[145., 122., 130., 95., 115.],
 #       [80., 63., 85., 48., 78.],
  #      [121., 107., 93., 69., 95.],
   #     [118., 83., 116., 80., 105.],
    #    [97., 75., 120., 80., 111.]]

cellTable, timeList = organizeData(timeData)

# intialize the objective function V (to be minimized) to a high value
# longest case time multiplied by number of cases >= possible function value

value = timeList[-1].time * len(cellTable[0])

# initialize null alloc
# the index of alloc corresponds to a lawyer
# the value at that index is the case assigned to that lawyer

allocation = [-1]*len(timeData)

# start algo at index = 0

(allocation, value) = optimalAlloc(cellTable, timeList, allocation, value, 0)

print('Case number by lawyer index: ')
print(allocation)
print('Objective Function Value: ')
print(value)
