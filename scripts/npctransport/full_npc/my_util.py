from __future__ import print_function

def range_inclusive(first, last, step=1):
        return range(first,last+step,step)

if __name__ == "__main__":
        assert(range_inclusive(2, 0, -1) == [2, 1, 0])
        assert(range_inclusive(0, 4)== [ 0, 1, 2, 3, 4 ])
        assert(range_inclusive(0, 40, 10)== [ 0, 10, 20, 30, 40 ])