import os
import sys
import shutil

try:
  os.environ['gPLUTO_DIR']
except KeyError:
  print('gPLUTO_DIR not defined. Setting it to the Current Directory')
  pluto_directory = os.getcwd()
  pass
else:
  pluto_directory = os.environ['gPLUTO_DIR']

sys.path.append(pluto_directory + '/Tools/Python/')
import menu
import configure
from make_problem import MakeProblem


def PlutoInterface(pluto_dir, string_mod, do_auto_update = False):
  work_dir = os.getcwd()
  interface_optval = ''
  interface_opts = ['Setup problem', 'Change makefile',
                    'Auto-update','Save Setup','Quit']

  if do_auto_update:
    MakeProblem(work_dir, pluto_dir, 1, 1)
    interface_optval = interface_opts[-1]   # set to "Quit" so it'll skip next loop

  while interface_optval != interface_opts[-1]:
    if (string_mod == ""): string_mod = "[None]"  
    menu.SetTitle( "Python setup (March 2025)","Working dir:     "+work_dir
                  +"\nPLUTO dir:       "+pluto_dir
                  +"\nEnabled Modules: "+string_mod)

    interface_optval = menu.Browse(interface_opts)

    if interface_optval == interface_opts[0]:   #  Setup problem
      if not os.path.exists(work_dir+'/init.cpp'):
        shutil.copy(pluto_dir+'/Src/Templates/init.cpp',work_dir+'/init.cpp')

      if not os.path.exists(work_dir+'/pluto.ini'):
        shutil.copy(pluto_dir+'/Src/Templates/pluto.ini',work_dir+'/pluto.ini')

      MakeProblem(work_dir, pluto_dir, 0, 1)

    if interface_optval == interface_opts[1]:  #  Change makefile
      MakeProblem(work_dir, pluto_dir, 1, 0)

    if interface_optval == interface_opts[2]:  #  Auto-update
      menu.Prompt('Press Enter to Update '+work_dir)
      MakeProblem(work_dir, pluto_dir, 1, 1)
      menu.Print ("Configuration up to date",sleep=0.75)
      break

    if interface_optval == interface_opts[3]: #  Save Setup
      # ------------------------------------
      #  get the path basename and make
      #  a tar archive
      # ------------------------------------

      menu.RestoreScreen()
      os.chdir(work_dir)
      default_name = os.path.basename(work_dir)
      os.system("clear")
      setup_name = input(" > Archive name ("+default_name+"): ")
      complete_backup = 0
      scrh       = input(" > Backup source files ? (n): ")
      if (setup_name == ''):
        setup_name = default_name

      if (scrh == "y"):
        complete_backup = 1

      setup_name = setup_name+".tar"
      ifail = os.system ("tar cvf "+setup_name+" *.cpp *.ini *.hpp")

      # -------------------------------
      #  back up source tree
      # -------------------------------

      if (complete_backup):
        ifail = ifail or os.system("tar rvf "+setup_name+" --directory="+pluto_dir+"/   Src/ ")
        ifail = ifail or os.system("tar rvf "+setup_name+" --directory="+pluto_dir+"/   Config/")
        ifail = ifail or os.system("tar rvf "+setup_name+" makefile")

      ifail = ifail or os.system ("gzip -f "+setup_name)
      if (ifail == 0):
        print (" > "+setup_name+".gz successfully created\n")
      else:
        print (" ! Error creating "+setup_name+".gz\n")

      sys.exit()

  if (menu.CursesIsActive()): menu.RestoreScreen()
  print("\n> Done.")
  sys.exit()

if __name__ == "__main__":   # starts from here
  auto_update = 0
  print("\n> Checking system architecture\n")
  configure.check(pluto_directory, 1)
  
  string_mod = ""            # String description on which module is enabled
  for x in sys.argv[1:]:     # check argument list

    if (x == "--auto-update"):
      auto_update = 1

    elif (x == "--get-arch"):
      sys.exit(1)
      break

    elif (x == "--help" or x == "-help"):
      print ("Usage: python $gPLUTO_DIR/setup.py [options]\n")
      print ("Here [options] can be:\n")
      print (" --auto-update     Run the python script in background.")
      print (" --no-curses       Disable ncurses library and use a")
      print ("                   simpler text-based menu.")
      print (" --with-fargo      Enable the FARGO-MHD module.")
      print (" --with-ho         Enable support for high-order methods.")
      sys.exit(1)

    elif (x == "--no-curses"):
      print("")

    elif (x == "--with-fargo"):
      print ("Enabling support for FARGO scheme")

    elif (x == "--with-poisson_solver"):
      string_mod += "[POISSON_SOLVER]" 
      print ("Enabling support for Poisson Solver")

    elif (x == "--with-ho"):
      print ("Enabling support for 4th-order")
      string_mod += "[HIGH_ORDER]" 

    elif (x == "--with-sb"):
      print ("Enabling support for shearing box module")
      string_mod += "[SHEARING_BOX]" 
      if '--with-fd' in sys.argv:
        print('! Incompatible modules, ',x,' +  --with-fd')
        sys.exit(1)

    else:
      print ("! Unrecognized option '",x,"'")
      sys.exit(1)

  print('\n> Loading PLUTO Interface...')

  if auto_update == 1:
    PlutoInterface(pluto_directory,string_mod, do_auto_update=True)
  else:
    PlutoInterface(pluto_directory,string_mod)
