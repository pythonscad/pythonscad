from pythonscad import *
import ast
import inspect

import sys
sys.path.insert(0, "/home/gsohler/tmp")

from PyQt6 import sip , QtWidgets, QtCore, QtGui


def parse_arguments(cls, args_text):
    """Parst einen bestehenden Aufruf-String (Inhalt zwischen den Klammern)
    generisch anhand der Parameter, die cls.__new__ tatsaechlich akzeptiert.

    Rueckgabe: (bekannt, unbekannt)
      bekannt   -- dict {parametername: wert} fuer alles, was zur Signatur passt
      unbekannt -- dict {name: wert} fuer Keyword-Argumente, die im Aufruf
                   standen, aber nicht Teil der Signatur sind (z.B. weil die
                   Klasse sich seither geaendert hat) -- werden beim
                   Zurueckschreiben mit uebernommen, statt stillschweigend
                   verworfen zu werden.
    """
    sig = inspect.signature(cls.__new__)
    param_names = [
        name for name, p in sig.parameters.items()
        if name != "cls"
        and p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)
    ]

    known = {}
    unknown = {}
    if not args_text.strip():
        return known, unknown

    try:
        tree = ast.parse(f"f({args_text})", mode="eval")
    except SyntaxError:
        return known, unknown

    call = tree.body

    for i, arg_node in enumerate(call.args):
        if i < len(param_names):
            try:
                known[param_names[i]] = ast.literal_eval(arg_node)
            except (ValueError, SyntaxError):
                pass

    for kw in call.keywords:
        try:
            value = ast.literal_eval(kw.value)
        except (ValueError, SyntaxError):
            continue
        if kw.arg in param_names:
            known[kw.arg] = value
        else:
            unknown[kw.arg] = value

    return known, unknown

HIT_RADIUS = 10  # Bildschirm-Pixel-Toleranz zum Anklicken (zoom-unabhaengig)
SCALE = 20        # Pixel pro Welteinheit bei zoom=1.0

class PolygonCanvas(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.points = []  # Welt-Koordinaten, unabhaengig von Zoom/Pan
        self.setFixedSize(500, 500)
        self.pan = QtCore.QPointF(250, 250)   # Ursprung startet in der Mitte

        self.drag_index = None
        self.zoom = 1.0
        self.panning = False
        self.pan_start_mouse = None
        self.pan_start_pan = None

        self.on_change = None

    def _to_world(self, screen_pos):
        return QtCore.QPointF(
            (screen_pos.x() - self.pan.x()) / self.zoom,
            -(screen_pos.y() - self.pan.y()) / self.zoom,   # Y invertiert
        )

    def _to_screen(self, world_xy):
        return QtCore.QPointF(
            world_xy[0] * self.zoom + self.pan.x(),
            -world_xy[1] * self.zoom + self.pan.y(),        # Y invertiert
        )

    def world_points(self):
        """Welt-Pixel -> Einheiten (SCALE=20 -> 20 = 1). Keine Hoehen-Abhaengigkeit mehr."""
        return [(px / SCALE, py / SCALE) for px, py in self.points]


    @staticmethod
    def _dist(a, b):
        return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5

    def _hit_test(self, x, y):
        radius_world = HIT_RADIUS / self.zoom
        for i, (px, py) in enumerate(self.points):
            if self._dist((px, py), (x, y)) <= radius_world:
                return i
        return None

    def _insert_at_best_edge(self, x, y):
        n = len(self.points)
        if n < 2:
            self.points.append((x, y))
            return
        best_idx, best_cost = None, None
        for i in range(n):
            p1 = self.points[i]
            p2 = self.points[(i + 1) % n]
            cost = self._dist(p1, (x, y)) + self._dist((x, y), p2) - self._dist(p1, p2)
            if best_cost is None or cost < best_cost:
                best_cost, best_idx = cost, i + 1
        self.points.insert(best_idx, (x, y))

    # --- Maus-Events ---
    def wheelEvent(self, event):
        old_world = self._to_world(event.position())
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        self.zoom = max(0.1, min(20.0, self.zoom * factor))
        new_screen = self._to_screen((old_world.x(), old_world.y()))
        self.pan += event.position() - new_screen  # Punkt unter dem Cursor bleibt fix
        self.update()

    def mousePressEvent(self, event):
        screen_pos = event.position()
        world_pos = self._to_world(screen_pos)
        hit = self._hit_test(world_pos.x(), world_pos.y())

        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            if hit is not None:
                self.drag_index = hit
            else:
                self._insert_at_best_edge(world_pos.x(), world_pos.y())
        elif event.button() == QtCore.Qt.MouseButton.RightButton:
            if hit is not None:
                del self.points[hit]
            else:
                self.panning = True
                self.pan_start_mouse = screen_pos
                self.pan_start_pan = QtCore.QPointF(self.pan)

        self.update()
        if self.on_change:
            self.on_change()

    def mouseMoveEvent(self, event):
        screen_pos = event.position()

        if self.drag_index is not None:
            world_pos = self._to_world(screen_pos)
            self.points[self.drag_index] = (world_pos.x(), world_pos.y())
            self.update()
            if self.on_change:
                self.on_change()
        elif self.panning:
            self.pan = self.pan_start_pan + (screen_pos - self.pan_start_mouse)
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self.drag_index = None
        elif event.button() == QtCore.Qt.MouseButton.RightButton:
            self.panning = False

    # --- Zeichnen ---
    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QtGui.QColor("white"))

        self._draw_axes(painter)

        screen_points = [self._to_screen(p) for p in self.points]

        if len(screen_points) > 1:
            pen_line = QtGui.QPen(QtGui.QColor("blue"))
            pen_line.setWidth(2)
            painter.setPen(pen_line)
            painter.drawPolyline(QtGui.QPolygonF(screen_points))
            painter.drawLine(screen_points[-1], screen_points[0])  # durchgezogene Schlusslinie

        for i, sp in enumerate(screen_points):
            color = "orange" if i == self.drag_index else "red"
            pen_point = QtGui.QPen(QtGui.QColor(color))
            pen_point.setWidth(8)
            painter.setPen(pen_point)
            painter.drawPoint(sp)

        painter.setPen(QtGui.QColor("black"))
        for i, sp in enumerate(screen_points):
            painter.drawText(int(sp.x()) + 8, int(sp.y()) - 8, str(i))

    def _draw_axes(self, painter):
        origin = self._to_screen((0, 0))
        w, h = self.width(), self.height()

        pen_x = QtGui.QPen(QtGui.QColor("red"))
        pen_x.setWidth(1)
        painter.setPen(pen_x)
        painter.drawLine(QtCore.QPointF(0, origin.y()), QtCore.QPointF(w, origin.y()))

        pen_y = QtGui.QPen(QtGui.QColor("blue"))
        pen_y.setWidth(1)
        painter.setPen(pen_y)
        painter.drawLine(QtCore.QPointF(origin.x(), 0), QtCore.QPointF(origin.x(), h))

        step_screen = SCALE * self.zoom
        if step_screen > 4:
            painter.setPen(QtGui.QColor("gray"))
            x = origin.x() % step_screen
            while x < w:
                painter.drawLine(QtCore.QPointF(x, origin.y() - 4), QtCore.QPointF(x, origin.y() + 4))
                x += step_screen
            y = origin.y() % step_screen
            while y < h:
                painter.drawLine(QtCore.QPointF(origin.x() - 4, y), QtCore.QPointF(origin.x() + 4, y))
                y += step_screen



_native_polygon = polygon
class polygon:
    @staticmethod
    def get_calltip():
        return "polygon(points=[[x1,y1], [x2,y2], ...])"

    @staticmethod
    def on_editor_trigger(pos):
        mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
        existing_args = editor_get_call_args(pos)
        known, unknown = parse_arguments(polygon, existing_args)

        dialog = QtWidgets.QDialog(mw)
        dialog.setWindowTitle("polygon")
        layout = QtWidgets.QVBoxLayout(dialog)

        hint = QtWidgets.QLabel(
            "Linksklick: Punkt hinzufuegen (rastet ein)   |   Rechtsklick: letzten Punkt entfernen"
        )
        layout.addWidget(hint)

        canvas = PolygonCanvas(dialog)
        world_pts = known.get("points", [])
        canvas.points = [(x * SCALE, y * SCALE) for x, y in world_pts]
        layout.addWidget(canvas)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        ok_button = buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        ok_button.setEnabled(len(canvas.points) >= 3)
        canvas.on_change = lambda: ok_button.setEnabled(len(canvas.points) >= 3)

        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            pts = canvas.world_points()
            points_str = ", ".join(f"[{x:.2f}, {y:.2f}]" for x, y in pts)
            parts = [f"points=[{points_str}]"]
            parts += [f"{k}={v!r}" for k, v in unknown.items()]  # unbekannte Args erhalten
            editor_replace_call_args(pos, ", ".join(parts))

    def __new__(cls, points):
        return _native_polygon(points)
