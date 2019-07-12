import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="NTMpy",
    version="0.0.2",
    author="Lukas Alber & Valentino Scalera",
    author_email="lukas.alber@fysik.su.se",
    description="A solver for parabolic partial differential equations",
    long_description="With this solver we provide an object which can be used to solve parabolic partial differential equations, i.e. 1st order in time and 2nd order in space.In addition we provide a class to compute the energy deposit by a source via Transfer Matrix Method and a class for graphical output for visualization.",
    long_description_content_type="text/markdown",
    url="https://github.com/udcm-su/heat-diffusion-1D",
    packages=setuptools.find_packages(),
    classifiers=[
	'Intended Audience :: Science/Research',
	'Topic :: Scientific/Engineering :: Physics',
    'Natural Language :: English',
    'Programming Language :: Python :: 3',
	'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    ],
    keywords='mathematical modeling two temperature model differential equation',
    install_requires = [
    'numpy' ,  
    'matplotlib',
    'bspline',
    'tqdm'
            ],
)

