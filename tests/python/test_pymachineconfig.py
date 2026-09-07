"""Unit tests for ``pymachineconfig.MachineConfig``.

Exercises the JSON-config and colour/power/feed bit-packing logic. The
``machineconfig`` C API call made by ``write_settings`` is a stubbed mock, so
no real settings are pushed. ``MachineConfig(name=None)`` builds the default
in-memory config without touching any config file.
"""

from __future__ import annotations

import os
import types

import pytest


@pytest.fixture
def mc(overlay: types.FunctionType) -> object:
    """Create a fresh ``MachineConfig(name=None)`` with clean class-level state."""
    pm = overlay("pymachineconfig")
    pm.MachineConfig._config = {}
    pm.MachineConfig._working = {}
    return pm.MachineConfig(name=None)


def test_default_config_has_color_table_and_default_entry(mc: object) -> None:
    """A default ``MachineConfig`` includes a colour table with known values
    (L00 = black, L02 = red) and a ``"default"`` section in the config dict.
    """
    assert mc.get_property_value("L00", "color") == 0x000000FF  # type: ignore[attr-defined]
    assert mc.get_property_value("L02", "color") == 0xFF0000FF  # type: ignore[attr-defined]
    assert "default" in mc.dict()  # type: ignore[attr-defined]


def test_get_types_and_labels_by_type(mc: object) -> None:
    """``get_types()`` returns known type names; ``get_label_by_type()`` returns their labels."""
    assert mc.get_types() == {"ColorTable", "default"}  # type: ignore[attr-defined]
    labels = mc.get_label_by_type("ColorTable")  # type: ignore[attr-defined]
    assert {"L00", "L02", "T1", "T2"} <= labels


def test_color2str_formats_hex(mc: object) -> None:
    """``color2str(L02)`` returns the hex string ``"#FF0000FF"`` for red."""
    assert mc.color2str("L02") == "#FF0000FF"  # type: ignore[attr-defined]


def test_gen_color_power_only_bit_packing(mc: object) -> None:
    """Power and feed values are packed into a 32-bit colour integer:
    power=1 at bit position 22 yields ``0x00400000``.
    """
    assert mc.gen_color(power=1, feed=-1) == 0x00400000  # type: ignore[attr-defined]
    assert mc.gen_color2hex(power=1, feed=-1) == "0x00400000"  # type: ignore[attr-defined]
    assert mc.gen_color2str(power=1, feed=-1) == "#00400000"  # type: ignore[attr-defined]


def test_gen_color_rejects_overrange_power(mc: object) -> None:
    """Power values exceeding the bit-field range return zero instead of wrapping."""
    assert mc.gen_color(power=1001, feed=-1) == 0  # type: ignore[attr-defined]


def test_configfile_absolute_path_passthrough(mc: object) -> None:
    """An absolute path passed to ``configfile()`` is returned unmodified."""
    assert mc.configfile("/etc/PythonSCAD.json") == "/etc/PythonSCAD.json"  # type: ignore[attr-defined]


def test_configfile_relative_resolves_under_home(mc: object, monkeypatch: pytest.MonkeyPatch) -> None:
    """A relative path is resolved under ``~/.config/PythonSCAD/``, respecting
    ``XDG_CONFIG_HOME`` and falling back to ``$HOME/.config`` when unset.
    """
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.setenv("HOME", "/tmp/pytest-home")
    resolved = mc.configfile("cfg.json")  # type: ignore[attr-defined]
    assert resolved == os.path.join("/tmp/pytest-home", ".config", "PythonSCAD", "cfg.json")


def test_write_then_read_json_roundtrip(mc: object, tmp_path: pytest.TempPathFactory) -> None:
    """Writing a config payload to disk and reading it back yields the same data."""
    target = str(tmp_path / "cfg.json")  # type: ignore[operator]
    payload = {"default": {"type": "default", "property": {"machine": None}}}
    mc.write(payload, name=target)  # type: ignore[attr-defined]
    assert mc.read(target) == payload  # type: ignore[attr-defined]


def test_set_property_value_updates_config(mc: object) -> None:
    """``set_property_value`` modifies a property and ``get_property_value`` returns the new value."""
    mc.set_property_value("L00", "power", 42)  # type: ignore[attr-defined]
    assert mc.get_property_value("L00", "power") == 42  # type: ignore[attr-defined]
