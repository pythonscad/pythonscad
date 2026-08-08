import threading
import time
from unittest.mock import patch

from setuptools import Distribution
from setuptools.command.build_ext import build_ext

from setup import BuildExtWithLexYacc, get_link_libraries


class FakeCompiler:
    def __init__(self):
        self.active = 0
        self.compiled = []
        self.lock = threading.Lock()
        self.max_active = 0

    def compile(self, sources, **kwargs):
        raise AssertionError("the command should replace this method")

    def _setup_compile(
        self,
        output_dir,
        macros,
        include_dirs,
        sources,
        depends,
        extra_postargs,
    ):
        objects = [f"{source}.o" for source in sources]
        build = dict(zip(objects, ((source, ".cc") for source in sources)))
        return macros, objects, extra_postargs, [], build

    def _get_cc_args(self, pp_opts, debug, extra_preargs):
        return []

    def _compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
        with self.lock:
            self.active += 1
            self.max_active = max(self.max_active, self.active)
        time.sleep(0.02)
        with self.lock:
            self.compiled.append(src)
            self.active -= 1


class MinimalCompiler:
    def __init__(self):
        self.compiled = []

    def compile(self, sources, **kwargs):
        self.compiled.extend(sources)


def test_build_ext_compiles_sources_in_parallel():
    sources = [f"source-{index}.cc" for index in range(6)]
    compiler = FakeCompiler()
    command = BuildExtWithLexYacc(Distribution())
    command.compiler = compiler
    command.parallel = 3

    def compile_extension(_command):
        compiler.compile(sources)

    with patch.object(build_ext, "build_extensions", compile_extension):
        command.build_extensions()

    assert compiler.max_active == 3
    assert sorted(compiler.compiled) == sources


def test_build_ext_uses_default_compile_without_private_compiler_hooks():
    sources = ["source.cc"]
    compiler = MinimalCompiler()
    command = BuildExtWithLexYacc(Distribution())
    command.compiler = compiler
    command.parallel = 3

    def compile_extension(_command):
        compiler.compile(sources)

    with patch.object(build_ext, "build_extensions", compile_extension):
        command.build_extensions()

    assert compiler.compiled == sources


def test_link_libraries_omit_unavailable_boost_system():
    with (
        patch("setup.get_pkg_config_libraries", return_value=[]),
        patch("setup.pkg_config_flags", return_value=[]),
        patch("setup.find_library", return_value=None),
        patch("setup.IS_DARWIN", False),
        patch("setup.IS_WINDOWS", False),
    ):
        libraries = get_link_libraries()

    assert "boost_regex" in libraries
    assert "boost_program_options" in libraries
    assert "boost_system" not in libraries


def test_link_libraries_include_available_boost_system():
    with (
        patch("setup.get_pkg_config_libraries", return_value=[]),
        patch("setup.pkg_config_flags", return_value=[]),
        patch("setup.find_library", return_value="libboost_system.so"),
        patch("setup.IS_DARWIN", False),
        patch("setup.IS_WINDOWS", False),
    ):
        libraries = get_link_libraries()

    assert "boost_system" in libraries
