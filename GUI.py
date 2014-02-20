import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import pylab
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc 
from itertools import product


from sympy import Matrix, Poly, symbols, pprint, zeros,sympify
#from matplotlibwidget import MatplotlibWidget
from check import check_input
from ALGO4 import ALGO4
from ALGO11_2 import ALGO11
from ALGO24 import ALGO24
from ALGO21 import ALGO21

from PySide import QtCore, QtGui
import sys

on_example=False  # in order to know if we are going to calculate an example . It helps to fill the matrices
already_entered={'A':False,'B':False,'C':False,'D':False} # in order to be able to check your input afterwards ---list maybe?
Names=['E','A','B','C','D']
As=zeros(2,2)
Bs=zeros(2,2)
Cs=zeros(2,2)
Ds=zeros(2,2)

class MainWindow(QtGui.QDialog):
    def __init__(self,  parent=None):
        super(MainWindow, self).__init__(parent)
        #QtGui.QDialog.__init__(self, parent)
        self.resize(500,550)
        tabWidget = QtGui.QTabWidget()
        tabWidget.addTab(InputTab(), self.tr("Input"))
        tabWidget.addTab(OutputTab(), self.tr("Output"))
       
        QuitButton=  QtGui.QPushButton(self.tr("Quit"))
   
        QuitButton.clicked.connect( self, QtCore.SLOT("close()"))
        #self.pushButton.clicked.connect(self.on_button_clicked) 
        #self.connect(cancelButton, QtCore.SIGNAL("clicked()"), self, QtCore.SLOT("reject()"))

        buttonLayout = QtGui.QHBoxLayout()
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(QuitButton)
        #buttonLayout.addWidget(cancelButton)

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(tabWidget)
        mainLayout.addLayout(buttonLayout)
        self.setLayout(mainLayout)

        self.setWindowTitle(self.tr("polynomial2gss"))
    
             
class InputTab(QtGui.QWidget):
    def __init__(self,parent=None):
    
        QtGui.QWidget.__init__(self, parent)
        
        #Create Buttons
        ButtonsGroup = QtGui.QGroupBox(self.tr("Input of Data:"))
        MatrixA = QtGui.QPushButton(self.tr("Matrix A"))
        MatrixB = QtGui.QPushButton(self.tr("Matrix B"))
        MatrixC = QtGui.QPushButton(self.tr("Matrix C"))
        MatrixD = QtGui.QPushButton(self.tr("Matrix D"))
        
        RadioGroup=QtGui.QGroupBox(self.tr("Type of Model:"))
        self.radio1 = QtGui.QRadioButton("Wolovich Equivalent Model")
        self.radio2 = QtGui.QRadioButton("Bosgra && Van Den Weiden  Model")
        self.radio3 = QtGui.QRadioButton("Tan && Vandewall  Model ")
        self.radio4 = QtGui.QRadioButton("Tan && Vandewall  Model (second algorithm)")
        
        #Buttons Layout
        ButtonsLayout = QtGui.QVBoxLayout()
        ButtonsLayout =QtGui.QVBoxLayout()
        ButtonsLayout.addWidget(MatrixA)
        ButtonsLayout.addWidget(MatrixB)
        ButtonsLayout.addWidget(MatrixC)
        ButtonsLayout.addWidget(MatrixD)
        
        ButtonsGroup.setLayout(ButtonsLayout)
        
        #Radiobuttons Layout
        RadioLayout=QtGui.QVBoxLayout()
        RadioLayout.addWidget(self.radio1)
        RadioLayout.addWidget(self.radio2)
        RadioLayout.addWidget(self.radio3)
        RadioLayout.addWidget(self.radio4)
        RadioGroup.setLayout(RadioLayout)
        
        
        Run = QtGui.QPushButton(self.tr("Run"))
        example1 = QtGui.QPushButton(self.tr("Load example 1"))
        example2 = QtGui.QPushButton(self.tr("Load example 2"))
        example3 = QtGui.QPushButton(self.tr("Load example 3"))
        example1.setToolTip('After pressing the button check the input Matrices  to see the example \'s input')
        example2.setToolTip('After pressing the button check the input Matrices  to see the example \'s input')
        example3.setToolTip('After pressing the button check the input Matrices  to see the example \'s input')
        self.checkBox = QtGui.QCheckBox('Include tests')

        #Connect Buttons
        MatrixA.clicked.connect(self.Open_DialogA)
        MatrixB.clicked.connect(self.Open_DialogB)
        MatrixC.clicked.connect(self.Open_DialogC)
        MatrixD.clicked.connect(self.Open_DialogD)
        Run.clicked.connect(self.onRun)
        example1.clicked.connect(self.example1)
        example2.clicked.connect(self.example2)
        example3.clicked.connect(self.example3)
        
        """
        grid = QtGui.QGridLayout()
        #Create Button Layout
        grid.addWidget(MatrixA, 1,1)
        grid.addWidget(MatrixB, 2,1)
        grid.addWidget(MatrixC, 3,1)
        grid.addWidget(MatrixD, 4,1)
        grid.addWidget(self.radio1,1, 2)
        grid.addWidget(self.radio2, 2,2)
        grid.addWidget(self.radio3, 3,2)
        grid.addWidget(self.radio4, 4,2)
        """
        #progress_bar #http://srinikom.github.io/pyside-docs/PySide/QtGui/QProgressBar.html
        self.progressBar = QtGui.QProgressBar()
        self.progressBar.setRange(0,1)

        
        #ButtonsLayout.addWidget(Run)
        #ButtonsGroup.setLayout(ButtonsLayout)
        #ButtonsGroup.setLayout(grid)
        #ButtonsLayout2 = QtGui.QVBoxLayout()
        #ButtonsLayout2.addWidget(Run)
        #Create Latex Layout
        LatexGroup = QtGui.QGroupBox(self.tr("Input of Data"))
        LatexLayout=QtGui.QVBoxLayout()
        #LatexLayout.addStretch(1)
        string=r'$A(s)\beta(t) =B(s)u(t)$' '\n'  r'$y(t)=C(s)\beta(t) +D(s)u(t) $' 
        LatexLayout.addWidget(Draw(self,string))
        LatexGroup.setLayout(LatexLayout)

        #add a horizontal line to slit Run button from others
        #horizontalLine =QtGui.QFrame()
        #horizontalLine.setFrameStyle(QtGui.QFrame.HLine)
        #horizontalLine.frameWidth=3
        
        #verticalLine =QtGui.QFrame()
        #verticalLine.setFrameStyle(QtGui.QFrame.VLine)
        #verticalLine.frameWidth=3
        
        #horizontalLine.lineWidth=3
        #horizontalLine.midlineWidth=3
        #horizontalLine.setSizePolicy(QSizePolicy.Minimum,QSizePolicy.Expanding)
        
        #Add everything to the main Layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(LatexGroup)
        
        inputLayout = QtGui.QHBoxLayout()
        inputLayout.addWidget(ButtonsGroup)
        #inputLayout.addWidget(verticalLine)
        inputLayout.addWidget(RadioGroup)
        mainLayout.addLayout(inputLayout)
        #mainLayout.addWidget(horizontalLine)
        
        grid = QtGui.QGridLayout()
        #Create Button Layout
        grid.addWidget(example1, 1,1)
        grid.addWidget(example2, 2,1)
        grid.addWidget(example3, 3,1)
        #grid.addWidget(MatrixC, 3,1)
        #grid.addWidget(MatrixD, 4,1)
        #grid.addWidget(self.radio1,1, 2)
        #grid.addWidget(self.radio2, 2,2)
        grid.addWidget(self.checkBox, 4,1)
        grid.addWidget(Run, 4,2)
        grid.addWidget(self.progressBar,5,2)
        #mainLayout.addWidget(Run)
        #mainLayout.addWidget(example1)
        #mainLayout.addWidget(example2)
        #mainLayout.addWidget(self.checkBox )
        mainLayout.addLayout(grid)
        self.setLayout(mainLayout)
        
    def Open_DialogA(self):
        
        self.dialog = Input_Dialog('A')
        self.dialog.exec_()
    
    def Open_DialogB(self):
        
        self.dialog = Input_Dialog('B')
        self.dialog.exec_()
           
    def Open_DialogC(self):
        
        self.dialog = Input_Dialog('C')
        self.dialog.exec_()
     
    def Open_DialogD(self):
        
        self.dialog = Input_Dialog('D')
        self.dialog.exec_()

    def onRun(self):
        global OUTPUT,on_example
        
        on_example=False # in order to prevent erasing matrices with old ones
        self.progressBar.reset()
        if check_input(As,Bs,Cs,Ds)==False:
            msgBox = QtGui.QMessageBox()
            msgBox.setWindowTitle(self.tr("polynomial2gss"))
            msgBox.setText("<b> Your input doesn't meet the problem requirements <\b>")
            msgBox.exec_()
        else:        
            if self.checkBox.isChecked():
                do_test=True
            else:
                do_test=False
        

            
            self.progressBar.setValue(0.1)
            if self.radio1.isChecked():
                OUTPUT=ALGO4(As,Bs,Cs,Ds,do_test)
            elif self.radio2.isChecked():
                OUTPUT=ALGO11(As,Bs,Cs,Ds,do_test)
            elif self.radio3.isChecked():
                OUTPUT=ALGO21(As,Bs,Cs,Ds,do_test)
            elif self.radio4.isChecked():
                OUTPUT=ALGO24(As,Bs,Cs,Ds,do_test)
            self.progressBar.setValue(1.0)
            
   

        
        

    def example1(self):
        global As,Bs,Cs,Ds, on_example
        on_example=True
        s=symbols('s')
        As=Matrix(2,2,[s+1 ,s**3+2*s**2,s**2 +3*s+2,s**4+4*s**3 +4*s**2 +s+2]);
        Bs=Matrix([[s**2+1],[s**3+2*s**2+s+3]]);
        Cs=Matrix([[-s**2-3*s-1,-s**4 -4*s**3-4*s**2+1]]);
        Ds=Matrix([[s**3+2*s**2 +s +2]]);
        msgBox = QtGui.QMessageBox()                                
        msgBox.setWindowTitle(self.tr("polynomial2gss"))
        msgBox.setText("<b> Push the Matrix Buttons to see the input matrices<\b>")
        msgBox.exec_()
    
    def example2(self):
        global As,Bs,Cs,Ds, on_example
        on_example=True
        s=symbols('s')
        As=Matrix(2,2,[s+1 ,s,s**2,s**4 ]);
        Bs=Matrix([[s+1],[s**5+2*s**2+s+3]]); # to 5 htan 3
        Cs=Matrix([[-s**2-3*s-1,1]]);
        Ds=Matrix([[s**3 +s +2]]);
        msgBox = QtGui.QMessageBox()                                
        msgBox.setWindowTitle(self.tr("polynomial2gss"))
        msgBox.setText("<b> Push the Matrix Buttons to see the input matrices<\b>")
        msgBox.exec_()
        
    def example3(self):
        global As,Bs,Cs,Ds, on_example
        on_example=True
        s=symbols('s')
        As=Matrix(3,3,[s**4+s**3+1,s+3,s**3+2,s**2+2,s**3+2,s**4+3,s**4+2,s**2+1,s-2])
        Bs=Matrix(3,2,[s**4+2,s-2,s**3+1,s-3,s**4+s**2+s+1,s**3+1])
        Cs=Matrix(1,3,[s**4+2,s**3+s+1,s**2+s+2])
        Ds=Matrix(1,2,[s**3+2,s**2+s+4])
        msgBox = QtGui.QMessageBox()                                
        msgBox.setWindowTitle(self.tr("polynomial2gss"))
        msgBox.setText("<b> Push the Matrix Buttons to see the input matrices<\b>")
        msgBox.exec_()
       
       

def Draw(self,string):
    """
    Draw Formula
    """
     # Get window background color
    bg = self.palette().window().color()
    cl = (bg.redF(), bg.greenF(), bg.blueF())
    fig = Figure(figsize=(5,1.2), dpi=80,facecolor=cl, edgecolor=cl)
    #http://stackoverflow.com/questions/2660319/putting-newline-in-matplotlib-label-with-tex-in-python  comment number 2
     
    #fig.text(0.25,0.35,string, fontsize=30)
    fig.suptitle(string, fontsize=30)
    # generate the canvas to display the plot
    canvas = FigureCanvas(fig)
    
    return canvas  
      

class OutputTab(QtGui.QWidget):
    def __init__(self,parent=None):
        QtGui.QWidget.__init__(self, parent)
        string=r'$E\dot{x}(t)=Ax(t) +Bu(t)$' '\n'  r'$y(t)=Cx(t) + Dy(t) $' 
        #LatexFrame
        LatexGroup = QtGui.QGroupBox(self.tr("Output"))
        LatexLayout=QtGui.QVBoxLayout()
        LatexLayout.addWidget(Draw(self,string))
        LatexGroup.setLayout(LatexLayout)
        
        #Button
        #ResultsButton = QtGui.QPushButton(self.tr("Results"))
        #ResultsButton.clicked.connect(self.Results)
        
        ViewEButton= QtGui.QPushButton(self.tr("View E"))
        ViewAButton= QtGui.QPushButton(self.tr("View A"))
        ViewBButton= QtGui.QPushButton(self.tr("View B"))
        ViewCButton= QtGui.QPushButton(self.tr("View C"))
        ViewDButton= QtGui.QPushButton(self.tr("View D"))
        Validate= QtGui.QPushButton(self.tr("Validate The result"))
        """
        ViewEButton.clicked.connect(Results(self,0))
        ViewAButton.clicked.connect(Results(self,1))
        ViewBButton.clicked.connect(Results(self,2))
        ViewCButton.clicked.connect(Results(self,3))
        """
        ViewEButton.clicked.connect(self.ResultsE)
        ViewAButton.clicked.connect(self.ResultsA)
        ViewBButton.clicked.connect(self.ResultsB)
        ViewCButton.clicked.connect(self.ResultsC)
        ViewDButton.clicked.connect(self.ResultsD)
        Validate.clicked.connect(self.ValidateT)
        #ViewAButton.clicked.connect(Results(self,1))
        #ViewBButton.clicked.connect(self.Results(2))
        #ViewCButton.clicked.connect(self.Results(3))
        #Layout Settings
        mainLayout = QtGui.QVBoxLayout()
        
        mainLayout.addWidget(LatexGroup)
        mainLayout.addWidget(ViewEButton)
        mainLayout.addWidget(ViewAButton)
        mainLayout.addWidget(ViewBButton)
        mainLayout.addWidget(ViewCButton)
        mainLayout.addWidget(ViewDButton)
        mainLayout.addWidget(Validate)
        #mainLayout.addWidget(ResultsButton)
        #mainLayout.addStretch(1)
        self.setLayout(mainLayout)
    
    def ResultsE(self):
        self.dialog=Output_Dialog(Names[0],OUTPUT[0])
        self.dialog.exec_()
    
    def ResultsA(self):
        self.dialog=Output_Dialog(Names[1],OUTPUT[1])
        self.dialog.exec_()
    
    def ResultsB(self):
        self.dialog=Output_Dialog(Names[2],OUTPUT[2])
        self.dialog.exec_()
    
    def ResultsC(self):
        self.dialog=Output_Dialog(Names[3],OUTPUT[3])
        self.dialog.exec_()
        
    def ResultsD(self):
        self.dialog=Output_Dialog(Names[4],OUTPUT[4])
        self.dialog.exec_()

    
    def ValidateT(self):
        
            if OUTPUT[5]==True:
                text="<b>Yes<b/>"
            elif OUTPUT[5]==False:
                text="<b>No</b>"
            else:
                text="<b>Results aren't validated</b>"
            
            msgBox = QtGui.QMessageBox()
            msgBox.setWindowTitle(self.tr("polynomial2gss"))
            msgBox.setText("Are the results fullfil the Furhmann system-equivalence ? :" )
            msgBox.setInformativeText(text)
            msgBox.exec_()

def Results(self,Matrix_index):

    i=Matrix_index
    self.dialog=Output_Dialog(Names[i],OUTPUT[i],parent=self)
    self.dialog.exec_()        
    
class Output_Dialog(QtGui.QDialog):
    def __init__(self,Name_of_Matrix,Data_Matrix,parent=None):
        QtGui.QDialog.__init__(self)
        self.tableWidget = QtGui.QTableWidget()
        self.tableWidget.setColumnCount(Data_Matrix.cols)
        self.tableWidget.setRowCount(Data_Matrix.rows)
        
        for i in range(Data_Matrix.rows):
            for j in range(Data_Matrix.cols):
                item = QtGui.QTableWidgetItem(str(Data_Matrix[i,j]))
                self.tableWidget.setItem(i, j, item)
                
        ExportButton = QtGui.QPushButton(self.tr("Export"))
        ExportButton.clicked.connect(self.export)
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.tableWidget)
        mainLayout.addWidget(ExportButton)
        self.setLayout(mainLayout)
        string='Matrix:'+str(Name_of_Matrix)
        self.resize(700,300)
        self.setWindowTitle(self.tr(string))
    
    def export(self):
        #http://stackoverflow.com/questions/1230222/selected-rows-line-in-qtableview-copy-to-qclipboard
        text=""
        for i in range(self.tableWidget.rowCount()):
            text +='\n'
            for j in range(self.tableWidget.columnCount()):
                text +=self.tableWidget.item(i,j).text()
                text +='\t'
        
        clipboard = QtGui.QApplication.clipboard()
        clipboard.setText(text)
        msgBox = QtGui.QMessageBox()
        msgBox.setWindowTitle(self.tr("polynomial2gss"))
        msgBox.setText("The results copied to system's clipboard !" )
        msgBox.exec_()

        
class Input_Dialog(QtGui.QDialog):
    def __init__(self, Name_of_Matrix, parent=None):
        QtGui.QDialog.__init__(self)
        
        self.Name_of_Matrix=str(Name_of_Matrix) #get the input matrix name
        #Create the table
        self.tableWidget = QtGui.QTableWidget()
        
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(2)

        self.fill()
        
        
        #Create Buttons 
        AddRowButton = QtGui.QPushButton(self.tr("Add Row"))
        RemoveRowButton=QtGui.QPushButton(self.tr("Remove Last Row"))
        AddColumnButton = QtGui.QPushButton(self.tr("Add Column"))
        RemoveColumnButton=QtGui.QPushButton(self.tr("Remove Last Column"))
        ApplyButton = QtGui.QPushButton(self.tr("Apply"))
        ApplyButton = QtGui.QPushButton(self.tr("Apply"))
        #Connect Buttons
        AddRowButton.clicked.connect(self.Add_Row)
        RemoveRowButton.clicked.connect(self.Remove_Row)
        AddColumnButton.clicked.connect(self.Add_Column)
        RemoveColumnButton.clicked.connect(self.Remove_Column)
        ApplyButton.clicked.connect(self.GET)
        """
        buttonLayout = QtGui.QHBoxLayout()
        buttonLayout.addStretch(1)
        """
        buttonLayout = QtGui.QGridLayout()
        buttonLayout.addWidget(AddRowButton,1,1)
        buttonLayout.addWidget(RemoveRowButton,2,1)
        buttonLayout.addWidget(AddColumnButton,1,2)
        buttonLayout.addWidget(RemoveColumnButton,2,2)
        buttonLayout.addWidget(ApplyButton,3,3)

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.tableWidget)
        mainLayout.addWidget(QtGui.QLabel('Use the symbol <b> s </b> as variable'))
        mainLayout.addLayout(buttonLayout)
        
        self.setLayout(mainLayout)
        string='Matrix: '+str(Name_of_Matrix)
        self.setWindowTitle(self.tr(string))
               
        
    def Add_Row(self):
        # http://zetcode.com/gui/pysidetutorial/eventsandsignals/
        row=self.tableWidget.rowCount()
        self.tableWidget.setRowCount(row +1)
    
    def Remove_Row(self):
        row=self.tableWidget.rowCount()
        self.tableWidget.setRowCount(row -1)
        
    def Add_Column(self):
        # http://zetcode.com/gui/pysidetutorial/eventsandsignals/
        col=self.tableWidget.columnCount()
        self.tableWidget.setColumnCount(col +1)
    
    def Remove_Column(self):
        # http://zetcode.com/gui/pysidetutorial/eventsandsignals/
        col=self.tableWidget.columnCount()
        self.tableWidget.setColumnCount(col -1)
    
    def GET(self):
        """
        Get elements from  Qtablewidget
        #http://stackoverflow.com/a/3771844/3026983
        """
        from itertools import product
        s=symbols('s')
        data=[sympify(self.tableWidget.item(i,j).text()) for i,j in  product(range(self.tableWidget.rowCount()), range(self.tableWidget.columnCount()))]
        
        if self.Name_of_Matrix=='A':
            global As,already_entered
            As=Matrix(self.tableWidget.rowCount(),self.tableWidget.columnCount(),data)
            already_entered['A']=True
            
        elif self.Name_of_Matrix=='B':
            global Bs ,already_entered
            Bs=Matrix(self.tableWidget.rowCount(),self.tableWidget.columnCount(),data) 
            already_entered['B']=True
            
        elif self.Name_of_Matrix=='C':
            global Cs,already_entered
            Cs=Matrix(self.tableWidget.rowCount(),self.tableWidget.columnCount(),data) 
            already_entered['C']=True
        else:
            global Ds ,already_entered
            Ds=Matrix(self.tableWidget.rowCount(),self.tableWidget.columnCount(),data) 
            already_entered['D']=True
        
        
        self.close()        #close input window
        
    def fill(self):
        if self.Name_of_Matrix=='A':
            M=As
        elif self.Name_of_Matrix=='B':
            M=Bs
            
        elif self.Name_of_Matrix=='C':
            M=Cs 
        else:
            M=Ds

        self.tableWidget.setColumnCount(M.cols)
        self.tableWidget.setRowCount(M.rows)
        for i,j in product(range(M.rows),range(M.cols)):
            item = QtGui.QTableWidgetItem(str(M[i,j]).replace('**','^'))  #sympy uses ** for exponents
            self.tableWidget.setItem(i, j, item)
        


        
        
