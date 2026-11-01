"""GUI smoke tests (offscreen Qt). They document the current GUI behaviour:
the main window opens, resizes, and runs a reconstruction in a QThread."""
import os
import time

import numpy as np
import pytest

pytest.importorskip('pytestqt')
qt = pytest.importorskip('pytomo.qt')
from pytomo.qt import QTimer, QMessageBox  # noqa: E402

from conftest import default_inputs  # noqa: E402


@pytest.fixture(scope='module')
def window(qapp, d3d_server):
    # One window per module: creating a second MainWindow in the same QApplication
    # currently segfaults (Qt object lifetime issue in GUI.py, see report).
    from conftest import _get_d3d_sxr, D3D_SHOT
    from pytomo import GUI
    mp = pytest.MonkeyPatch()
    # never block on modal dialogs
    mp.setattr(QMessageBox, 'question', staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    mp.setattr(QMessageBox, 'warning', staticmethod(lambda *a, **k: QMessageBox.StandardButton.Ok))
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


def test_binding(window):
    print('Qt binding:', qt.API)
    assert qt.API == os.environ.get('PYTOMO_QT', 'PyQt6')


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


def test_cancel_reconstruction(window, qtbot):
    """STOP must interrupt the worker between blocks (was a no-op: wrong Qt import)."""
    import pytomo.main as M
    calls = []
    orig = M.interruption_requested
    def requested():
        calls.append(1)
        return orig()
    M.interruption_requested = requested
    try:
        window.thread.requestInterruption()  # before start: QThread resets it on start
        window.actionStart()
        window.thread.requestInterruption()
        qtbot.waitUntil(lambda: not window.thread.isRunning(), timeout=180_000)
    finally:
        M.interruption_requested = orig
    assert calls  # the worker polled for interruption


def test_toggle_all_controls(window, qtbot):
    """Every checkbox and combo box (and their slots) must work with the Qt binding in use."""
    from pytomo.qt import QCheckBox, QComboBox
    boxes = window.findChildren(QCheckBox)
    combos = window.findChildren(QComboBox)
    assert len(boxes) > 5 and len(combos) > 3
    for b in boxes:
        state = b.isChecked()
        b.setChecked(not state)
        b.setChecked(state)
    for c in combos:
        i = c.currentIndex()
        for j in range(c.count()):
            c.setCurrentIndex(j)
        c.setCurrentIndex(i)
    qtbot.wait(50)


def test_data_setting_window(window, qtbot):
    """The data selection / preprocessing window (data_preprocessing.py) opens and closes."""
    from pytomo.data_preprocessing import DataSettingWindow
    window.showSetting()
    qtbot.wait(200)
    wins = [w for w in window.findChildren(DataSettingWindow)]
    assert wins and wins[-1].isVisible()
    for w in wins:
        w.close()
    qtbot.wait(50)
