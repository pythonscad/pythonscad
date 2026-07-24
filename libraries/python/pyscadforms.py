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

HIT_RADIUS = 10
DRAG_THRESHOLD = 4  # Pixel Bewegung, ab der ein Klick als Rubber-Band-Drag zaehlt
SCALE = 20


class PolygonCanvas(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.points = []          # Welt-Koordinaten
        self.selected = set()     # Indizes der ausgewaehlten Punkte
        self.setFixedSize(500, 500)
        self.pan = QtCore.QPointF(250, 250)

        self.zoom = 1.0
        self.panning = False
        self.pan_start_mouse = None
        self.pan_start_pan = None

        self.drag_active = False
        self.drag_start_world = None
        self.drag_orig_points = None  # {index: (x,y)} beim Drag-Start

        self.press_screen_pos = None
        self.press_was_hit = None
        self.press_shift = False
        self.rubber_band_active = False
        self.rubber_band_rect = None  # QRectF in Bildschirmkoordinaten

        self.on_change = None          # Punkte haben sich geaendert (Werte)
        self.on_selection_change = None  # Auswahl hat sich geaendert

    # --- Koordinatentransformation ---
    def _to_world(self, screen_pos):
        return QtCore.QPointF(
            (screen_pos.x() - self.pan.x()) / self.zoom,
            -(screen_pos.y() - self.pan.y()) / self.zoom,
        )

    def _to_screen(self, world_xy):
        return QtCore.QPointF(
            world_xy[0] * self.zoom + self.pan.x(),
            -world_xy[1] * self.zoom + self.pan.y(),
        )

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
            return len(self.points) - 1
        best_idx, best_cost = None, None
        for i in range(n):
            p1 = self.points[i]
            p2 = self.points[(i + 1) % n]
            cost = self._dist(p1, (x, y)) + self._dist((x, y), p2) - self._dist(p1, p2)
            if best_cost is None or cost < best_cost:
                best_cost, best_idx = cost, i + 1
        self.points.insert(best_idx, (x, y))
        return best_idx

    # --- Aktionen, auch von aussen (Tabelle/Buttons) aufrufbar ---
    def delete_selected(self):
        if not self.selected:
            return
        self.points = [p for i, p in enumerate(self.points) if i not in self.selected]
        self.selected = set()
        self._notify_change()
        self._notify_selection()

    def add_point_default(self):
        if self.points:
            last = self.points[-1]
            new_pt = (last[0] + 1.0, last[1] + 1.0)
        else:
            new_pt = (0.0, 0.0)
        idx = self._insert_at_best_edge(*new_pt)
        self.selected = {idx}
        self._notify_change()
        self._notify_selection()

    def set_point_value(self, index, x, y):
        if 0 <= index < len(self.points):
            self.points[index] = (x, y)
            self._notify_change()

    def set_selection(self, indices):
        self.selected = set(indices)
        self.update()
        self._notify_selection(from_table=True)

    def _notify_change(self):
        self.update()
        if self.on_change:
            self.on_change()

    def _notify_selection(self, from_table=False):
        self.update()
        if self.on_selection_change and not from_table:
            self.on_selection_change()

    # --- Maus-Events ---
    def wheelEvent(self, event):
        old_world = self._to_world(event.position())
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        self.zoom = max(0.1, min(20.0, self.zoom * factor))
        new_screen = self._to_screen((old_world.x(), old_world.y()))
        self.pan += event.position() - new_screen
        self.update()

    def mousePressEvent(self, event):
        screen_pos = event.position()
        world_pos = self._to_world(screen_pos)
        hit = self._hit_test(world_pos.x(), world_pos.y())
        shift = bool(event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier)

        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            if hit is not None:
                if shift:
                    if hit in self.selected:
                        self.selected.discard(hit)
                    else:
                        self.selected.add(hit)
                    self._notify_selection()
                else:
                    if hit not in self.selected:
                        self.selected = {hit}
                        self._notify_selection()
                    self.drag_active = True
                    self.drag_start_world = world_pos
                    self.drag_orig_points = {i: self.points[i] for i in self.selected}
            else:
                # Weder Punkt getroffen -> entweder spaeter ein Klick (= Punkt
                # einfuegen) oder ein Rubber-Band-Drag. Entscheidet sich in
                # mouseMoveEvent/mouseReleaseEvent anhand der Bewegung.
                self.press_screen_pos = screen_pos
                self.press_shift = shift
                if not shift:
                    self.selected = set()
                    self._notify_selection()
        elif event.button() == QtCore.Qt.MouseButton.RightButton:
            if hit is not None:
                self.points.pop(hit)
                self.selected.discard(hit)
                self.selected = {i - 1 if i > hit else i for i in self.selected}
                self._notify_change()
                self._notify_selection()
            else:
                self.panning = True
                self.pan_start_mouse = screen_pos
                self.pan_start_pan = QtCore.QPointF(self.pan)

        self.update()
        if self.on_change:
            self.on_change()

    def mouseMoveEvent(self, event):
        screen_pos = event.position()

        if self.drag_active:
            world_pos = self._to_world(screen_pos)
            dx = world_pos.x() - self.drag_start_world.x()
            dy = world_pos.y() - self.drag_start_world.y()
            for i, (ox, oy) in self.drag_orig_points.items():
                self.points[i] = (ox + dx, oy + dy)
            self._notify_change()
        elif self.panning:
            self.pan = self.pan_start_pan + (screen_pos - self.pan_start_mouse)
            self.update()
        elif self.press_screen_pos is not None:
            delta = screen_pos - self.press_screen_pos
            if not self.rubber_band_active and (abs(delta.x()) > DRAG_THRESHOLD or abs(delta.y()) > DRAG_THRESHOLD):
                self.rubber_band_active = True
            if self.rubber_band_active:
                self.rubber_band_rect = QtCore.QRectF(self.press_screen_pos, screen_pos).normalized()
                self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            if self.rubber_band_active and self.rubber_band_rect is not None:
                hit_indices = {
                    i for i, p in enumerate(self.points)
                    if self.rubber_band_rect.contains(self._to_screen(p))
                }
                if self.press_shift:
                    self.selected |= hit_indices
                else:
                    self.selected = hit_indices
                self._notify_selection()
            elif self.press_screen_pos is not None:
                # Reiner Klick ohne Bewegung -> Punkt einfuegen
                world_pos = self._to_world(self.press_screen_pos)
                idx = self._insert_at_best_edge(world_pos.x(), world_pos.y())
                self.selected = {idx}
                self._notify_change()
                self._notify_selection()

            self.drag_active = False
            self.press_screen_pos = None
            self.rubber_band_active = False
            self.rubber_band_rect = None
            self.update()
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
            painter.drawLine(screen_points[-1], screen_points[0])

        for i, sp in enumerate(screen_points):
            color = "orange" if i in self.selected else "red"
            pen_point = QtGui.QPen(QtGui.QColor(color))
            pen_point.setWidth(8)
            painter.setPen(pen_point)
            painter.drawPoint(sp)

        painter.setPen(QtGui.QColor("black"))
        for i, (sp, (wx, wy)) in enumerate(zip(screen_points, self.points)):
            painter.drawText(int(sp.x()) + 8, int(sp.y()) - 8, f"({wx:.2f}, {wy:.2f})")

        if self.rubber_band_active and self.rubber_band_rect is not None:
            pen_rb = QtGui.QPen(QtGui.QColor(0, 120, 255))
            pen_rb.setStyle(QtCore.Qt.PenStyle.DashLine)
            painter.setPen(pen_rb)
            painter.setBrush(QtGui.QColor(0, 120, 255, 40))
            painter.drawRect(self.rubber_band_rect)

    def _draw_axes(self, painter):
        origin = self._to_screen((0, 0))
        w, h = self.width(), self.height()

        pen_x = QtGui.QPen(QtGui.QColor("red")); pen_x.setWidth(1)
        painter.setPen(pen_x)
        painter.drawLine(QtCore.QPointF(0, origin.y()), QtCore.QPointF(w, origin.y()))

        pen_y = QtGui.QPen(QtGui.QColor("blue")); pen_y.setWidth(1)
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

    def world_points(self):
        return [(px / SCALE, py / SCALE) for px, py in self.points]


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
        outer = QtWidgets.QVBoxLayout(dialog)

        hint = QtWidgets.QLabel(
            "Klick auf leere Flaeche: Punkt einfuegen (an guenstigster Stelle)\n"
            "Klick+Ziehen auf Punkt: verschieben (auch mehrere ausgewaehlte)\n"
            "Shift+Klick: Punkt zur Auswahl hinzufuegen/entfernen\n"
            "Ziehen auf leerer Flaeche: Mehrfachauswahl per Rahmen\n"
            "Rechtsklick auf Punkt: einzeln loeschen   |   Rechtsklick leer: Ansicht verschieben"
        )
        outer.addWidget(hint)

        body = QtWidgets.QHBoxLayout()
        outer.addLayout(body)

        canvas = PolygonCanvas(dialog)
        world_pts = known.get("points", [])
        canvas.points = [(x * SCALE, y * SCALE) for x, y in world_pts]
        body.addWidget(canvas)

        side = QtWidgets.QVBoxLayout()
        body.addLayout(side)

        table = QtWidgets.QTableWidget(0, 2)
        table.setHorizontalHeaderLabels(["X", "Y"])
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        side.addWidget(table)

        btn_add = QtWidgets.QPushButton("Punkt hinzufuegen")
        btn_delete = QtWidgets.QPushButton("Auswahl loeschen")
        side.addWidget(btn_add)
        side.addWidget(btn_delete)

        syncing = {"active": False}  # verhindert Rueckkopplung zwischen Tabelle <-> Canvas

        def refresh_table():
            syncing["active"] = True
            table.setRowCount(len(canvas.points))
            for i, (x, y) in enumerate(canvas.world_points()):
                table.setItem(i, 0, QtWidgets.QTableWidgetItem(f"{x:.3f}"))
                table.setItem(i, 1, QtWidgets.QTableWidgetItem(f"{y:.3f}"))
            table.blockSignals(False)
            syncing["active"] = False

        def refresh_table_selection():
            syncing["active"] = True
            table.clearSelection()
            for i in canvas.selected:
                table.selectRow(i)
            syncing["active"] = False

        def on_cell_changed(row, col):
            if syncing["active"]:
                return
            try:
                x = float(table.item(row, 0).text()) * SCALE
                y = float(table.item(row, 1).text()) * SCALE
                canvas.set_point_value(row, x, y)
            except (ValueError, AttributeError):
                pass

        def on_table_selection_changed():
            if syncing["active"]:
                return
            rows = {idx.row() for idx in table.selectionModel().selectedRows()}
            canvas.set_selection(rows)

        canvas.on_change = refresh_table
        canvas.on_selection_change = refresh_table_selection
        table.cellChanged.connect(on_cell_changed)
        table.itemSelectionChanged.connect(on_table_selection_changed)
        btn_add.clicked.connect(canvas.add_point_default)
        btn_delete.clicked.connect(canvas.delete_selected)

        refresh_table()

        delete_shortcut = QtGui.QShortcut(QtGui.QKeySequence(QtCore.Qt.Key.Key_Delete), dialog)
        delete_shortcut.activated.connect(canvas.delete_selected)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        outer.addWidget(buttons)

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            pts = canvas.world_points()
            points_str = ", ".join(f"[{x:.2f}, {y:.2f}]" for x, y in pts)
            parts = [f"points=[{points_str}]"]
            parts += [f"{k}={v!r}" for k, v in unknown.items()]
            editor_replace_call_args(pos, ", ".join(parts))

    def __new__(cls, points):
        return _native_polygon(points)

class CubePreview(QtWidgets.QWidget):
    # Asymmetrische Projektionswinkel statt symmetrischer Isometrie (30/30),
    # damit ein Wuerfel nicht wie sechs gleiche Dreiecke wirkt.
    ANGLE_X_DEG = 20
    ANGLE_Y_DEG = 40

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedSize(220, 220)
        self.size = (1.0, 1.0, 1.0)

    def set_size(self, x, y, z):
        self.size = (max(x, 0.001), max(y, 0.001), max(z, 0.001))
        self.update()

    @classmethod
    def _iso_raw(cls, dx, dy, dz):
        import math
        ax = math.radians(cls.ANGLE_X_DEG)
        ay = math.radians(cls.ANGLE_Y_DEG)
        ix = dx * math.cos(ax) - dy * math.cos(ay)
        iy = -dz - dx * math.sin(ax) - dy * math.sin(ay)
        return ix, iy

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QtGui.QColor("white"))

        x, y, z = self.size
        max_dim = max(x, y, z)
        target = 100
        scale = target / max_dim
        sx, sy, sz = x * scale, y * scale, z * scale

        raw = {}
        for key in ("000", "100", "010", "001", "110", "101", "011", "111"):
            dx = sx if key[0] == "1" else 0
            dy = sy if key[1] == "1" else 0
            dz = sz if key[2] == "1" else 0
            raw[key] = self._iso_raw(dx, dy, dz)

        xs = [p[0] for p in raw.values()]
        ys = [p[1] for p in raw.values()]
        bbox_w = max(xs) - min(xs) or 1
        bbox_h = max(ys) - min(ys) or 1
        mid_x = (max(xs) + min(xs)) / 2
        mid_y = (max(ys) + min(ys)) / 2

        pad = 30
        avail_w = self.width() - 2 * pad
        avail_h = self.height() - 2 * pad
        fit_scale = min(avail_w / bbox_w, avail_h / bbox_h, 1.0)

        cx, cy = self.width() / 2, self.height() / 2
        corners = {
            key: QtCore.QPointF(cx + (rx - mid_x) * fit_scale, cy + (ry - mid_y) * fit_scale)
            for key, (rx, ry) in raw.items()
        }

        # Kanten je nach Richtung (X/Y/Z) farblich zugeordnet
        edges_x = [("000", "100"), ("010", "110"), ("001", "101"), ("011", "111")]
        edges_y = [("000", "010"), ("100", "110"), ("001", "011"), ("101", "111")]
        edges_z = [("000", "001"), ("100", "101"), ("010", "011"), ("110", "111")]

        for edges, color in ((edges_x, "red"), (edges_y, "darkgreen"), (edges_z, "blue")):
            pen = QtGui.QPen(QtGui.QColor(color))
            pen.setWidth(2)
            painter.setPen(pen)
            for a, b in edges:
                painter.drawLine(corners[a], corners[b])

        painter.setPen(QtGui.QColor("red"))
        p = corners["100"]
        painter.drawText(int(p.x()) - 10, int(p.y()) + 18, f"x={x:.2f}")

        painter.setPen(QtGui.QColor("darkgreen"))
        p = corners["010"]
        painter.drawText(int(p.x()) - 45, int(p.y()) + 4, f"y={y:.2f}")

        painter.setPen(QtGui.QColor("blue"))
        p = corners["001"]
        painter.drawText(int(p.x()) + 6, int(p.y()) - 6, f"z={z:.2f}")

_native_cube = cube

class cube:
    @staticmethod
    def get_calltip():
        return "cube(size=[1,1,1], center=False)"

    def __new__(cls, size=1, center=False):
        return _native_cube(size, center)

    @staticmethod
    def on_editor_trigger(pos):
        mw = sip.wrapinstance(mainwindow_ptr(), QtWidgets.QMainWindow)
        existing_args = editor_get_call_args(pos)
        known, unknown = parse_arguments(cube, existing_args)

        dialog = QtWidgets.QDialog(mw)
        dialog.setWindowTitle("cube")
        outer = QtWidgets.QHBoxLayout(dialog)

        preview = CubePreview(dialog)
        outer.addWidget(preview)

        form_col = QtWidgets.QVBoxLayout()
        outer.addLayout(form_col)

        size_val = known.get("size", 1)
        if not isinstance(size_val, (list, tuple)):
            size_val = [size_val, size_val, size_val]

        form = QtWidgets.QFormLayout()
        form_col.addLayout(form)

        x_box = QtWidgets.QDoubleSpinBox(); x_box.setRange(0.01, 10000); x_box.setValue(size_val[0])
        y_box = QtWidgets.QDoubleSpinBox(); y_box.setRange(0.01, 10000); y_box.setValue(size_val[1])
        z_box = QtWidgets.QDoubleSpinBox(); z_box.setRange(0.01, 10000); z_box.setValue(size_val[2])
        form.addRow("size x", x_box)
        form.addRow("size y", y_box)
        form.addRow("size z", z_box)

        uniform_box = QtWidgets.QCheckBox("Einheitliche Groesse (X/Y/Z gekoppelt)")
        form_col.addWidget(uniform_box)

        center_box = QtWidgets.QCheckBox()
        form.addRow("center", center_box)
        center_box.setChecked(bool(known.get("center", False)))

        preview.set_size(x_box.value(), y_box.value(), z_box.value())

        _syncing = {"active": False}

        def on_value_changed(changed_box):
            if _syncing["active"]:
                return
            if uniform_box.isChecked():
                _syncing["active"] = True
                value = changed_box.value()
                x_box.setValue(value)
                y_box.setValue(value)
                z_box.setValue(value)
                _syncing["active"] = False
            preview.set_size(x_box.value(), y_box.value(), z_box.value())

        x_box.valueChanged.connect(lambda _: on_value_changed(x_box))
        y_box.valueChanged.connect(lambda _: on_value_changed(y_box))
        z_box.valueChanged.connect(lambda _: on_value_changed(z_box))

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        form_col.addWidget(buttons)

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            size_str = f"[{x_box.value():.3f}, {y_box.value():.3f}, {z_box.value():.3f}]"
            parts = [f"size={size_str}", f"center={center_box.isChecked()}"]
            parts += [f"{k}={v!r}" for k, v in unknown.items()]
            editor_replace_call_args(pos, ", ".join(parts))
