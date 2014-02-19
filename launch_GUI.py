import sys
"""
from PySide import QtGui, QtCore
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import pylab
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc 
"""
from GUI import *
"""
def main():
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = MainWindow()
    ui.show()  # ?? command eixei mia
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()
"""


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    tabdialog = MainWindow()
    sys.exit(tabdialog.exec_())
