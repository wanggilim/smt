from setuptools import setup, find_packages
from sphinx.setup_command import BuildDoc

cmdclass = {'build_sphinx':BuildDoc}
name = 'sof-dossier'
lname = "SOFIA Dossier Tool"
version = "1.0.1"
release = version


setup(
    name=name,
    version=version,
    author="Michael S. Gordon",
    author_email="msgordon.astro@gmail.com",
    description="Generate SOFIA dossiers from mission and flight IDs",
    packages=find_packages(),
    package_data={
        '':['*.tex','*.cfg']
    },
    install_requires=['astropy','beautifulsoup4','peewee',
                      'numpy','pandas','mechanicalsoup','scipy',
                      'urlpath','matplotlib','regions','aplpy',
                      'requests','astroquery','shapely','pylatexenc'],
    cmdclass=cmdclass,
    command_options={
        'build_sphinx': {
            'project': ('setup.py', lname),
            'version': ('setup.py', version),
            'release': ('setup.py', release),
            'source_dir': ('setup.py', 'doc')}},
    entry_points={
        'console_scripts': [
            'sof-dossier = dossier.dossier:main',
            'sof-planner = dossier.planner:main',
        ]
    },
    zip_safe=False
)
