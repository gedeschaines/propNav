# -*- coding: utf-8 -*-

# pylint: disable=trailing-whitespace,bad-whitespace,invalid-name
# pylint: disable=anomalous-backslash-in-string,bad-continuation
# pylint: disable=multiple-statements,redefined-outer-name,global-statement

"""
FILE:  pqueLib.py
DATE:  06 DEC 2023
AUTH:  G. E. Deschaines
DESC:  Data structures and methods for a priority queue implemented
       with an ordered binary heap maintained in an array. A decent
       explanation of priority queues can be found at this link:
       
         http://algs4.cs.princeton.edu/24pq/

       Methods presented herein were derived from those presented in
       a textbook on programming and data structures using Pascal.    
REFS:
    
  [1] This Python script was refactored from pquelib.c available at:
      https://github.com/gedeschaines/threeD/blob/master/src/pquelib.c


Disclaimer:  See DISCLAIMER file.
    
"""

class HeapElement:
    
    def __init__(self, key=0, info=0):
        self.key  = key
        self.info = info


class PriorityQueue:
    
    def __init__(self, maxElements):
        self._maxElements = maxElements
        self._bottom = 0
        self._Elements = [HeapElement(0, k) for k in range(maxElements+1)]

    # Private methods.
    
    def _reheapUp(self, bottom):
        
        heapOk       = False
        currentIndex = bottom
        parentIndex  = currentIndex // 2
        while (currentIndex > 1) and (not heapOk):
            if ( self._Elements[parentIndex].key >= 
                 self._Elements[currentIndex].key ):
                heapOk = True
            else:
                tempElement                  = self._Elements[parentIndex]
                self._Elements[parentIndex]  = self._Elements[currentIndex]
                self._Elements[currentIndex] = tempElement
                currentIndex                 = parentIndex
                parentIndex                  = parentIndex // 2

    def _reheapDown(self, root, bottom):

        root2   = root*2
        root2p1 = root2 + 1
        heapOk  = False
        while (root2 <= bottom) and (not heapOk):
            if root2 == bottom: 
                maxChild = root2
            else:
                if self._Elements[root2].key > self._Elements[root2p1].key:
                    maxChild = root2
                else:
                    maxChild = root2p1
            if self._Elements[root].key < self._Elements[maxChild].key:
                tempElement              = self._Elements[root]
                self._Elements[root]     = self._Elements[maxChild]
                self._Elements[maxChild] = tempElement
                root                     = maxChild
                root2                    = root*2
                root2p1                  = root2 + 1
            else:
                heapOk = True

    # Public methods.
    
    def clear(self):
        self._bottom = 0
        
    def isEmpty(self):
        return self._bottom == 0

    def isFull(self):
        return self._bottom == self._maxElements
    
    def priorityEnq(self, newElement):
        self._bottom                 = self._bottom + 1
        self._Elements[self._bottom] = newElement
        self._reheapUp(self._bottom)

    def priorityDeq(self):
        firstElement      = self._Elements[1]
        self._Elements[1] = self._Elements[self._bottom]
        self._bottom      = self._bottom - 1
        self._reheapDown(1, self._bottom)
        return firstElement
