"""Qt binding used by the GUI: PyQt6, or PyQt5 when that is already loaded
(pyTomo windows embedded in a PyQt5 application, e.g. pyspecview) or PyQt6 is missing.

Code using this module must use fully scoped enums (Qt.AlignmentFlag.AlignCenter),
which both bindings accept. Set PYTOMO_QT=PyQt5 to force PyQt5.
"""
import os
import sys

API = os.environ.get('PYTOMO_QT')
if API is None:
    API = 'PyQt5' if 'PyQt5.QtCore' in sys.modules else 'PyQt6'

if API == 'PyQt6':
    try:
        from PyQt6 import QtCore, QtGui, QtWidgets
        from PyQt6.QtGui import QAction
    except ImportError:
        API = 'PyQt5'
if API == 'PyQt5':
    from PyQt5 import QtCore, QtGui, QtWidgets
    from PyQt5.QtWidgets import QAction

# matplotlib selects the binding through QT_API (backend_qtagg supports both)
os.environ.setdefault('QT_API', API.lower())

Signal = QtCore.pyqtSignal
Slot = QtCore.pyqtSlot
Qt = QtCore.Qt
QObject = QtCore.QObject
QRect = QtCore.QRect
QThread = QtCore.QThread
QTimer = QtCore.QTimer
QIcon = QtGui.QIcon
QPalette = QtGui.QPalette
QPixmap = QtGui.QPixmap
QApplication = QtWidgets.QApplication
QCheckBox = QtWidgets.QCheckBox
QComboBox = QtWidgets.QComboBox
QDoubleSpinBox = QtWidgets.QDoubleSpinBox
QFileDialog = QtWidgets.QFileDialog
QGridLayout = QtWidgets.QGridLayout
QGroupBox = QtWidgets.QGroupBox
QHBoxLayout = QtWidgets.QHBoxLayout
QLabel = QtWidgets.QLabel
QLineEdit = QtWidgets.QLineEdit
QMainWindow = QtWidgets.QMainWindow
QMessageBox = QtWidgets.QMessageBox
QProgressBar = QtWidgets.QProgressBar
QPushButton = QtWidgets.QPushButton
QSizePolicy = QtWidgets.QSizePolicy
QSlider = QtWidgets.QSlider
QSpacerItem = QtWidgets.QSpacerItem
QSpinBox = QtWidgets.QSpinBox
QStyleFactory = QtWidgets.QStyleFactory
QTabWidget = QtWidgets.QTabWidget
QToolButton = QtWidgets.QToolButton
QVBoxLayout = QtWidgets.QVBoxLayout
QWidget = QtWidgets.QWidget


def center_on_screen(widget):
    """Move a top-level widget to the centre of its screen."""
    screen = widget.screen() if hasattr(widget, 'screen') else QtWidgets.QApplication.primaryScreen()
    frame = widget.frameGeometry()
    frame.moveCenter(screen.availableGeometry().center())
    widget.move(frame.topLeft())
