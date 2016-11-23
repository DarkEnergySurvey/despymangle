import distutils
from distutils.core import setup
import glob

bin_files = glob.glob("bin/*")

# The main call
setup(name='despymangle',
      version ='2.2.4',
      license = "GPL",
      description = "DES mangle framework",
      author = "Aurelien Benoit-Levy, Molly Swanson, Michelle Gower",
      author_email = "ucapab2@ucl.ac.uk",
      packages = ['despymangle'],
      package_dir = {'': 'python'},
      scripts = bin_files,
      data_files=[('ups',['ups/despymangle.table'])]
      )
