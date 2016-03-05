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

    # data table and lists by lawyer created

    # now create sorted list of times in ascending order for additive algorithm

    # sortedList contains Cell objects

    sortedList = []

    for i in range(len(data)):
        for j in range(len(data[0])):
            sortedList.append(table[i][j])
    
    sortedList = sorted(sortedList, key=lambda Cell: Cell.time)

    return table, sortedList


def updateIncumbent(table, alloc):

    # calculate the updated allocation's obj func value

    newValue = 0

    for i in range(len(alloc)):
        newValue += table[i][alloc[i]].time

    return newValue


def updateAlloc(table):

    # update alloc from cells in table

    alloc = [-1]*len(timeData)
    for i in range(len(table)):
        for j in range(len(table[0])):
            if table[i][j].assigned == 1:
                alloc[i] = j

    return alloc

# recursively determine the optimal allocation of cases
# using Balas additive algorithm, assume root node cannot be solution

def optimalAlloc(table, sortList, index, prevAlloc, prevIncumbent):

    # update return values

    alloc = updateAlloc(table)
    incumbent = updateIncumbent(table, alloc)

    # <index> is the index of the binary variable in the sorted obj func
    
    # check constraints (one case per lawyer)

    # node is impossible if 2 cases have been assigned to a lawyer

    # node is infeasible if 0 cases have been assigned to a lawyer
    # (all lawyers must be given a case s.t. 2 cases are not assigned
    # to a lawyer, because size(lawyers)=size(cases))

    # if a node is feasible, save the allocation and obj func value,
    # then recurse until every branch is pruned

    # in the lawyer problem, the feasible nodes will be in seperate branches,
    # so we can return if a feasible node is found

    # this simplification is generic to binary knapsack problems
    # i.e. one-to-one mapping of element sets

    # list of number of cases per lawyer, indexed by lawyer

    sumList = []

    for i in range(len(table)):
        tempSum = 0
        for j in range(len(table[0])):
            tempSum += table[i][j].assigned

        if tempSum > 1:  # lawyer i assigned to an impossible number of cases
            # exit recursion with the previous node's feasible allocation
            # (prune node)

            return (prevAlloc, prevIncumbent)

        sumList.append(tempSum)

    # Exit recursion if this is the last variable, so as to not recurse again
    if len(sortList)-1 == index:
        return (alloc, incumbent)

    for i in sumList:
        # Check for infeasible constraint

        # if node is infeasible, we must branch

        # if there are no cases assigned for some lawyer,
        # this node is infeasible
        if i == 0:
            # Depth First Recursion

            # options A and B are the branches off
            # of the current node in the binary tree
            # the two opetions are chosen by choosing to assign
            # either the current case (A) or the next case (B) to lawyer i

            # assign the current cell in the sorted list
            table[sortList[index].lawyer][sortList[index].case].assigned = 1

            # recurse option A, evaluate the next index
            (allocA, incumbentA) = optimalAlloc(table, sortList, index+1, alloc, incumbent)

            # assign the next cell in the sorted list
            table[sortList[index].lawyer][sortList[index].case].assigned = 0
            table[sortList[index+1].lawyer][sortList[index+1].case].assigned = 1

            # recurse option B, evaluate the next index
            (allocB, incumbentB) = optimalAlloc(table, sortList, index+1, alloc, incumbent)

            # save the alloc with the better obj func value
            # in this lawyer problem, "better" is less time
            if incumbentB < incumbentA:
                return (allocB, incumbentB)
            return (allocA, incumbentA)

    # satisfied constraints

    return (alloc, incumbent)



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


# start algo at index = 0

startIndex = 0

# initialize the previous return values

prevAlloc = updateAlloc(cellTable)

prevIncumbent = updateIncumbent(cellTable, prevAlloc)

# run algorithm

(allocation, value) = optimalAlloc(cellTable, timeList, startIndex, prevAlloc, prevIncumbent)

print('Case number by lawyer index: ')
print(allocation)
print('Objective Function Value: ')
print(value)

