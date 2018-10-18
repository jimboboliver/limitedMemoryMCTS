from PyQt5 import QtWidgets, QtGui, QtCore
from visualiser import Ui_MainWindow
import sys
from subprocess import Popen, PIPE
import subprocess
from time import sleep


class PhotoViewer(QtWidgets.QGraphicsView):
    photoClicked = QtCore.pyqtSignal(QtCore.QPoint)

    def __init__(self, parent):
        super(PhotoViewer, self).__init__(parent)
        self._zoom = 0
        self._empty = True
        self._scene = QtWidgets.QGraphicsScene(self)
        self._photo = QtWidgets.QGraphicsPixmapItem()
        self._scene.addItem(self._photo)
        self.setScene(self._scene)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor(30, 30, 30)))
        self.setFrameShape(QtWidgets.QFrame.NoFrame)

    def hasPhoto(self):
        return not self._empty

    def fitInView(self, scale=True):
        rect = QtCore.QRectF(self._photo.pixmap().rect())
        if not rect.isNull():
            self.setSceneRect(rect)
            if self.hasPhoto():
                unity = self.transform().mapRect(QtCore.QRectF(0, 0, 1, 1))
                self.scale(1 / unity.width(), 1 / unity.height())
                viewrect = self.viewport().rect()
                scenerect = self.transform().mapRect(rect)
                factor = min(viewrect.width() / scenerect.width(),
                             viewrect.height() / scenerect.height())
                self.scale(factor, factor)
            self._zoom = 0

    def setPhoto(self, pixmap=None):
        self._zoom = 0
        if pixmap and not pixmap.isNull():
            self._empty = False
            self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
            self._photo.setPixmap(pixmap)
        else:
            self._empty = True
            self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
            self._photo.setPixmap(QtGui.QPixmap())
        self.fitInView()

    def wheelEvent(self, event):
        if self.hasPhoto():
            if event.angleDelta().y() > 0:
                factor = 1.25
                self._zoom += 1
            else:
                factor = 0.8
                self._zoom -= 1
            if self._zoom > 0:
                self.scale(factor, factor)
            elif self._zoom == 0:
                self.fitInView()
            else:
                self._zoom = 0

    def toggleDragMode(self):
        if self.dragMode() == QtWidgets.QGraphicsView.ScrollHandDrag:
            self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
        elif not self._photo.pixmap().isNull():
            self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)

    def mousePressEvent(self, event):
        if self._photo.isUnderMouse():
            self.photoClicked.emit(QtCore.QPoint(event.pos()))
        super(PhotoViewer, self).mousePressEvent(event)


class window(QtWidgets.QMainWindow):
    iteration = -1
    max = -1
    states = []
    regens = []
    total_regens = []
    ref_regens = []
    total_ref_regens = []
    init = True

    def __init__(self):
        self.p = Popen('./visualmcts', stdin=PIPE, stdout=PIPE, encoding='utf8')

        super(window, self).__init__()

        self.ui = Ui_MainWindow()

        self.ui.setupUi(self)

        self.ui.nextButton.clicked.connect(self.next)
        self.ui.nextButton.setAutoDefault(True)
        self.ui.nextButton.setFocus()
        self.ui.previousButton.clicked.connect(self.previous)
        self.ui.previousButton.setAutoDefault(True)
        self.ui.playButton.clicked.connect(self.make_play)
        self.ui.playEdit.returnPressed.connect(self.make_play)

        self.imageDisplay = PhotoViewer(self)
        self.ui.verticalLayout.replaceWidget(self.ui.graphicsView, self.imageDisplay)
        self.ui.graphicsView.deleteLater()

        self.ui.playEdit.setEnabled(False)
        self.ui.playButton.setEnabled(False)
        self.ui.gameLabel.setText("***\n***\n***" + '\n\n')

    def reset(self):
        self.iteration = -1
        self.max = -1
        subprocess.call(['rm graph/graph*'], shell=True)
        self.states.clear()
        self.regens.clear()
        self.total_regens.clear()
        self.ref_regens.clear()
        self.total_ref_regens.clear()

    def update(self):
        self.ui.iterationNumber.setText(str(self.iteration + 1))

        if self.init:  # fudge, otherwise the first image is set to be tiny
            self.init = False
            self.imageDisplay.setPhoto(QtGui.QPixmap('graph/graph{}.jpg'.format(self.iteration)))
            sleep(0.001)
        self.imageDisplay.setPhoto(QtGui.QPixmap('graph/graph{}.jpg'.format(self.iteration)))

        self.ui.statesNumber.setText(str(self.states[self.iteration]))
        self.ui.regenNumber.setText(str(self.regens[self.iteration]))
        self.ui.regenTotalNumber.setText(str(self.total_regens[self.iteration]))
        self.ui.regenRefNumber.setText(str(self.ref_regens[self.iteration]))
        self.ui.regenRefTotalNumber.setText(str(self.total_ref_regens[self.iteration]))

    def next(self):
        self.iteration += 1
        if self.iteration > self.max:
            read = self.p.stdout.readline()
            cleaned = read
            if '*' in read or 'X' in read or 'O' in read:
                read = self.p.stdout.read(21)
                cleaned += read
            if 'Enter play' not in read and 'Game over' not in read:
                if self.iteration > self.max:
                    self.p.stdin.write("\n")
                    self.p.stdin.flush()
                    subprocess.check_call(["dot", 'graph/graph{}.dot'.format(self.iteration), "-Tjpg", "-o", 'graph/graph{}.jpg'.format(self.iteration)])

                    states_info = [int(x) for x in cleaned.split()]
                    self.states.append(states_info[0])
                    self.regens.append(states_info[1])
                    if len(self.total_regens) == 0:
                        self.total_regens.append(states_info[1])
                    else:
                        self.total_regens.append(self.total_regens[-1] + states_info[1])
                    self.ref_regens.append(states_info[2])
                    if len(self.total_ref_regens) == 0:
                        self.total_ref_regens.append(states_info[2])
                    else:
                        self.total_ref_regens.append(self.total_ref_regens[-1] + states_info[2])

            else:
                self.reset()
                if 'Game over' not in read:
                    self.ui.playButton.setEnabled(True)
                    self.ui.playEdit.setEnabled(True)
                    self.ui.playEdit.setFocus()
                else:
                    cleaned += self.p.stdout.readline()
                self.ui.previousButton.setEnabled(False)
                self.ui.nextButton.setEnabled(False)
                self.ui.gameLabel.setText(cleaned.strip())
                return

        self.update()

        if self.max < self.iteration:
            self.max = self.iteration

    def previous(self):
        if self.iteration > 0:
            self.iteration -= 1
            self.update()

    def make_play(self):
        if self.ui.playEdit.text():
            self.p.stdin.write(self.ui.playEdit.text() + '\n')
            self.p.stdin.flush()

            read = self.p.stdout.read(2)
            if read[1] == 'C':
                self.p.stdout.read(29)
                self.ui.playText.setText("Try again!")
            else:
                self.ui.playText.setText("")
                self.ui.playEdit.clear()
                self.ui.playEdit.setEnabled(False)
                self.ui.previousButton.setEnabled(True)
                self.ui.nextButton.setEnabled(True)
                self.ui.nextButton.setAutoDefault(True)
                self.ui.nextButton.setFocus()
                self.ui.playButton.setEnabled(False)
                read += self.p.stdout.read(12)
                self.ui.gameLabel.setText(read.strip() + '\n\n')


    def closeEvent(self, event):
        self.p.kill()  # Make sure to kill the mcts subprocess on exit
        print("Killed mcts subprocess")
        event.accept()


subprocess.call(['rm graph/graph*'], shell=True)

app = QtWidgets.QApplication([])

application = window()

application.show()

sys.exit(app.exec())
