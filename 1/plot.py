import pylab
import numpy
import sys
import matplotlib.pyplot as plt


x_number_values = []
y_number_values = []
x_border = []
y_border = []
data = []

folder = 'pics/'
grids = 'grids/'

# Считывание хода метода
def inputGrid(path):
    global x_number_values
    global y_number_values
    global data
    x_number_values = []
    y_number_values = []
    data = []
    with open(path, 'r') as f:
        for line in f: # read rest of lines
            data.append([ float(x) for x in line.split()])
    for i in range(len(data)):
        x_number_values.append(data[i][0])
        y_number_values.append(data[i][1])


# Считывание хода метода
def inputBorder(path):
    global x_border
    global y_border
    global data
    x_border = []
    y_border = []
    data = []
    with open(path, 'r') as f:
        for line in f: # read rest of lines
            data.append([ float(x) for x in line.split()])
    for i in range(len(data)):
        x_border.append(data[i][0])
        y_border.append(data[i][1])


# Отрисовка хода метода
def draw_line(path):
    #plt.plot(x_border, y_border, linewidth=1)
    plt.title(path, fontsize=19)
    plt.xlabel('X', fontsize=10)
    plt.ylabel('Y', fontsize=10)
    plt.tick_params(axis='both', labelsize=8)

    

def makeData():
    x = numpy.arange (1, 256, 0.1)
    y = numpy.arange (1, 64, 0.1)
    xgrid, ygrid = numpy.meshgrid(x, y)
    zgrid = xgrid + ygrid
    return xgrid, ygrid, zgrid


def drawGrid(grid, border):
    inputGrid(grid)
    inputBorder(border)
    #x, y, z = makeData()
    #cs = pylab.contour(x, y, z, 25)
    draw_line(grid)
    for i in range(len(x_number_values)):
        plt.scatter(x_number_values[i], y_number_values[i], s=10, marker='s', color='grey')
    for i in range(len(x_border)):
        plt.scatter(x_border[i], y_border[i], s=10, marker='s', color='black')
    plt.savefig(folder + grid[6:] + '.png')
    #plt.show()
    plt.clf()

    


if __name__ == '__main__':
    grid = str(sys.argv[1])
    border = str(sys.argv[2])
    drawGrid(grid, border)