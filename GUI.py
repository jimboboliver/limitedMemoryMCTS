from PyQt5 import QtWidgets, QtGui, QtCore
from visualiser import Ui_MainWindow
import sys
from subprocess import Popen, PIPE
import subprocess


class Label(QtWidgets.QLabel):
    def __init__(self, img):
        super(Label, self).__init__()
        self.setFrameStyle(QtWidgets.QFrame.StyledPanel)
        self.pixmap = QtGui.QPixmap(img)

    def paintEvent(self, event):
        size = self.size()
        painter = QtGui.QPainter(self)
        point = QtCore.QPoint(0,0)
        scaledPix = self.pixmap.scaled(size, QtCore.Qt.KeepAspectRatio, transformMode = QtCore.Qt.SmoothTransformation)
        # start painting the label from left upper corner
        point.setX((size.width() - scaledPix.width())/2)
        point.setY((size.height() - scaledPix.height())/2)
        painter.drawPixmap(point, scaledPix)

    def changePixmap(self, img):
        self.pixmap = QtGui.QPixmap(img)
        self.repaint()


class window(QtWidgets.QMainWindow):
    iteration = -1
    max = -1
    init = False

    def __init__(self):
        self.p = Popen('./mcts', stdin=PIPE, encoding='utf8')

        super(window, self).__init__()

        self.ui = Ui_MainWindow()

        self.ui.setupUi(self)

        self.ui.nextButton.clicked.connect(self.next)
        self.ui.previousButton.clicked.connect(self.previous)

    def next(self):
        self.iteration += 1
        if self.iteration > self.max:
            self.p.stdin.write("\n")
            self.p.stdin.flush()
            subprocess.check_call(["dot", 'graph/graph{}.dot'.format(self.iteration), "-Tjpg", "-o", 'graph/graph{}.jpg'.format(self.iteration)])

        if not self.init:
            self.init = True
            self.imageLabel = Label('graph/graph{}.jpg'.format(self.iteration))
            self.ui.verticalLayout.replaceWidget(self.ui.imageLabel, self.imageLabel)

        self.imageLabel.changePixmap('graph/graph{}.jpg'.format(self.iteration))
        if self.max < self.iteration:
            self.max = self.iteration

    def previous(self):
        if self.iteration > 0:
            self.iteration -= 1
            self.imageLabel.changePixmap('graph/graph{}.jpg'.format(self.iteration))


subprocess.call(['rm graph/graph*'], shell=True)

app = QtWidgets.QApplication([])

application = window()

application.show()

sys.exit(app.exec())
