# -*- coding: cp1251 -*-
import pylab
import numpy
import sys
import matplotlib.pyplot as plt


x_inside = []
x_border = []
data = []

folder = 'pics/'
grids = 'grids/'

# Считывание хода метода
def inputGrid(path):
	global x_inside
	global data
	x_inside = []
	data = []
	with open(path, 'r') as f:
		for line in f: # read rest of lines
			data.append([ float(x) for x in line.split()])
	for i in range(len(data)):
		x_inside.append(data[i][0])


# Считывание хода метода
def inputBorder(path):
	global x_border
	global data
	x_border = []
	data = []
	with open(path, 'r') as f:
		for line in f: # read rest of lines
			data.append([ float(x) for x in line.split()])
	for i in range(len(data)):
		x_border.append(data[i][0])


# Отрисовка хода метода
def draw_line(title):
	#plt.plot(x_border, y_border, linewidth=1)
	plt.title(title, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	plt.tick_params(axis='both', labelsize=8)



def makeData():
	x = numpy.arange (1, 10, 0.1)
	y = numpy.arange (1, 10, 0.1)
	xgrid, ygrid = numpy.meshgrid(x, y)
	zgrid = xgrid + ygrid
	return xgrid, ygrid, zgrid


def drawGrid(grid, border):
	name = 'title'
	if grid[6:9] == 'Non':
		name = grid[6:16]
	else:
		name = grid[6:13]

	name +=' '
	inputGrid(grid)
	inputBorder(border)
	#x, y, z = makeData()
	#cs = pylab.contour(x, y, z, 25)
	draw_line(name + str(len(x_inside)+len(x_border))+' nodes')
	for i in range(len(x_inside)):
		plt.scatter(x_inside[i], 0, s=4, marker='s', color='grey')
	for i in range(len(x_border)):
		plt.scatter(x_border[i], 0, s=4, marker='s', color='black')
	plt.savefig(folder + grid[6:] + '.png')
	#plt.show()
	plt.clf()



if __name__ == '__main__':
	grid = str(sys.argv[1])
	border = str(sys.argv[2])
	drawGrid(grid, border)
