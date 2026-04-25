from setuptools import setup, find_packages

setup(
    name="pythonscad-stubs",
    version="0.2.0",
    description="Type stubs for PythonSCAD (_openscad / openscad / pythonscad)",
    packages=find_packages(),
    package_data={
        "_openscad": ["py.typed", "__init__.pyi"],
        "openscad": ["py.typed", "__init__.pyi"],
        "pythonscad": ["py.typed", "__init__.pyi"],
    },
    zip_safe=False,  # Required for type hints to work
    python_requires=">=3.11",
)
