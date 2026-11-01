"""GUI smoke tests (offscreen Qt). They document the current GUI behaviour:
the main window opens, resizes, and runs a reconstruction in a QThread."""
import time

import numpy as np
import pytest

pytest.importorskip('pytestqt')
QtWidgets = pytest.importorskip('PyQt5.QtWidgets')
from PyQt5.QtCore import QTimer  # noqa: E402

from conftest import default_inputs  # noqa: E402


@pytest.fixture(scope='module')
def window(qapp, d3d_server):
    # One window per module: creating a second MainWindow in the same QApplication
    # currently segfaults (Qt object lifetime issue in GUI.py, see report).
    from conftest import _get_d3d_sxr, D3D_SHOT
    from pytomo import GUI
    mp = pytest.MonkeyPatch()
    # never block on modal dialogs
    mp.setattr(GUI.QMessageBox, 'question', staticmethod(lambda *a, **k: GUI.QMessageBox.Yes))
    mp.setattr(GUI.QMessageBox, 'warning', staticmethod(lambda *a, **k: GUI.QMessageBox.Ok))
    _, tok = _get_d3d_sxr(D3D_SHOT)
    inputs = default_inputs(solver=3, regularization=1, plot_all=False, postprocessing=False,
                            screen_width=1920, screen_height=1080)
    inputs['tokamak_tmp'] = tok
    win = GUI.MainWindow(inputs, tok)
    win.show()
    yield win
    win.close()
    win.deleteLater()  # let Qt delete the C++ objects before Python's GC does (avoids exit segfault)
    qapp.processEvents()
    mp.undo()


def test_main_window_opens(window):
    tabs = [window.main_tab.tabText(i) for i in range(window.main_tab.count())]
    assert 'Main' in tabs
    assert window.windowTitle() == 'Tomography'


@pytest.mark.parametrize('size', [(800, 500), (1600, 1000)])
def test_resize_does_not_crash(window, qtbot, size):
    window.resize(*size)
    qtbot.wait(50)
    assert window.isVisible()


def test_reconstruction_in_thread(window, qtbot):
    """Runs START, measures GUI event-loop latency while the worker thread computes."""
    latencies = []
    last = [time.perf_counter()]

    def tick():
        now = time.perf_counter()
        latencies.append(now - last[0])
        last[0] = now

    timer = QTimer()
    timer.timeout.connect(tick)
    timer.start(20)
    window.actionStart()
    assert window.thread.isRunning()
    qtbot.waitUntil(lambda: not window.thread.isRunning(), timeout=180_000)
    qtbot.wait(200)
    timer.stop()
    pix = window.Picture_recons_1.pixmap()
    assert pix is not None and not pix.isNull()
    # documented, not enforced: worst event-loop stall during the computation
    print('max GUI stall %.0f ms' % (1e3 * max(latencies)))
