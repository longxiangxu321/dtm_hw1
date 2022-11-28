
import math
import numpy

class Raster:
    """A raster"""
    
    def __init__(self, cellsize, bbox):
        self.cellsize = cellsize
        self.bbox = bbox
        self.width = math.ceil((self.bbox[2] - self.bbox[0]) / cellsize)
        self.height = math.ceil((self.bbox[3] - self.bbox[1]) / cellsize)
        self.data = numpy.zeros(shape=(self.height, self.width))

    # @property
    # def cellsize(self):
    #     return self.cellsize
    # @cellsize.setter
    # def cellsize(self, value):
    #     self.value = value

    # @property
    # def width(self):
    #     return self.width

    # @property
    # def height(self):
    #     return self.height

    def __getitem__(self, indices):
        return self.data[indices]

    def __setitem__(self, indices, value):
        self.data[indices] = value

    def save(self, output_file):
        fout = open(output_file, 'w')
        fout.write('NCOLS %d\n' % self.data.shape[1])
        fout.write('NROWS %d\n' % self.data.shape[0])
        fout.write('XLLCORNER %d\n' % self.bbox[0])
        fout.write('YLLCORNER %d\n' % self.bbox[1])
        fout.write('CELLSIZE %d\n' % self.cellsize)
        fout.write('NODATA_VALUE -99999\n')
        for row in self:
            for each in row:
                fout.write(str(each) + ' ')
            fout.write('\n')
        fout.close()
        print("Raster succesfully saved as {}".format(output_file))

