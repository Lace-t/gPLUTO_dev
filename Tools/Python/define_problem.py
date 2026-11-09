import os
import sys
import menu as menu
import pluto_files_IO as pfIO

class DefineProblem(object):
  def __init__(self, work_dir, pluto_dir, auto_update):
    """Defines the problem for the PLUTO code.

    This class essentially creates the definitions.hpp file
    based on the inputs given by the user from the PlutoSetup menu.
    In case the header file is already present, this class will
    read those default values to re-create the file.

    **Inputs**:
      1. work_dir = Path to PLUTO code working directory
      2. pluto_dir = Path to PLUTO source directory
      3. auto-update = Boolean that indicates auto-update of definitions.hpp file.

    **Output**:
      It generates a definitions.hpp file.

    """
    #Some class attributes that will be used in various class methods.
    self.work_dir         = work_dir
    self.ini_fname        = self.work_dir + '/init.cpp'
    self.pluto_dir        = pluto_dir
    self.krome_dir        = pluto_dir + '/Src/Cooling/KROME/krome/'
    self.auto_update      = auto_update
    self.def_fname        = self.work_dir + '/definitions.hpp'
    self.additional_files = []
    self.header_files     = []
    self.additional_flags = []
    self.pluto_path       = []
    self.def_file_list    = []

    self.ChkCompatiblity() #To Replace Old Keywords with New Ones

    # defining the PLUTO entries and its default values in lists.
    self.entries = ['PHYSICS', 'DIMENSIONS', 'GEOMETRY','MAPPED_GRID',
                    'BODY_FORCE', 'COOLING', 'RECONSTRUCTION',
                    'TIME_STEPPING','NTRACER', 'PARTICLES','USER_DEF_PARAMETERS']
    self.default = ['HD', '1', 'CARTESIAN', 'NO',
                    'NO', 'NO','LINEAR',
                    'RK2', '0', 'NO','0']

    # Creating a dictionary of flags that are invoked by giving arguments.
    flag_keys = ['FULL', 'WITH-FD', 'WITH-SB', 'WITH-FARGO',
                 'WITH-HO','WITH-POISSON_SOLVER']
    #self.flag_dict = {key: False for key in flag_keys} DOESNT WORK WITH PYTHON 2.6
    self.flag_dict = {'WITH-CHOMBO':False, 'FULL':False, 'WITH-FD':False,
                      'WITH-SB':False, 'WITH-FARGO':False,
                      'WITH-HO':False, 'WITH-POISSON_SOLVER':False}

    for arg in sys.argv:
      if arg[2:].upper() in flag_keys:
        self.flag_dict[arg[2:].upper()] = True
      else:
        pass

    # Generates Full Option List.
    self.GenerateOptionsList()

    #Updates Options, default based on FLAGS.
    if True in self.flag_dict.values():
      self.AfterFlagLists()

    #Read the exsisting definition.hpp file or Browse the menu for Setting up problem.
    self.ReadOrBrowse(Ents=self.entries, Defs=self.default, Opts=self.options,
                      MenuTitle="Setup problem")

    #Process the PHYSICS Module.
    if self.default[self.entries.index('PHYSICS')] == 'HD':
      self.ProcessHDModule()
      self.ReadOrBrowse(Ents = self.entries_HD, Defs = self.default_HD,
                        Opts = self.options_HD, MenuTitle = "HD Menu")
      self.eos = self.default_HD[self.entries_HD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'RHD':
      self.ProcessRHDModule()
      self.ReadOrBrowse(Ents = self.entries_RHD, Defs = self.default_RHD,
                        Opts = self.options_RHD, MenuTitle = "RHD Menu")
      self.eos = self.default_RHD[self.entries_RHD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'MHD':
      self.ProcessMHDModule()
      self.ReadOrBrowse(Ents = self.entries_MHD, Defs = self.default_MHD,
                        Opts = self.options_MHD, MenuTitle = "MHD Menu")
      self.eos = self.default_MHD[self.entries_MHD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'RMHD':
      self.ProcessRMHDModule()
      self.ReadOrBrowse(Ents = self.entries_RMHD, Defs = self.default_RMHD,
                        Opts = self.options_RMHD, MenuTitle = "RMHD Menu")
      self.eos = self.default_RMHD[self.entries_RMHD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'ResRMHD':
      self.ProcessResRMHDModule()
      self.ReadOrBrowse(Ents = self.entries_ResRMHD, Defs = self.default_ResRMHD,
                        Opts = self.options_ResRMHD, MenuTitle = "ResRMHD Menu")
      self.eos = self.default_ResRMHD[self.entries_ResRMHD.index('EOS')]
    
    if self.default[self.entries.index('PHYSICS')] == 'GRMHD':
      self.ProcessGRMHDModule()
      self.ReadOrBrowse(Ents = self.entries_GRMHD, Defs = self.default_GRMHD,
                        Opts = self.options_GRMHD, MenuTitle = "GRMHD Menu")
      self.eos = self.default_GRMHD[self.entries_GRMHD.index('EOS')]
      
    # Process the KROME cooling module when required.
    if self.default[self.entries.index('COOLING')] == 'KROME':
      self.ProcessKROMEModule()
      self.ReadOrBrowse(Ents = self.entries_KROME, Defs = self.default_KROME,
                        Opts = self.options_KROME, MenuTitle = "KROME Menu")
      self.WriteKromeOpts()

    if self.default[self.entries.index('PARTICLES')] != 'NO':
      self.ProcessParticlesModule()
      self.ReadOrBrowse(Ents = self.entries_Particles, Defs = self.default_Particles,
                        Opts = self.options_Particles, MenuTitle = "PARTICLES Menu")

    #UserDef Parameters and Constants
    self.udef_params = []
    self.udef_const = []
    self.udef_const_vals = []
    self.nparam = int(self.default[self.entries.index('USER_DEF_PARAMETERS')])

    #self.nconst = int(self.default[self.entries.index('USER_DEF_CONSTANTS')])
    if self.nparam > 0:
      self.ProcessUserDefPara()

    self.ProcessUserDefConst()

    # Write a List def_file_list which will be written as the header file.
    self.WriteDefFileList()
    pf = pfIO.PlutoFiles(self.work_dir+'/definitions.hpp')
    pf.List2File(self.def_file_list)

  def ChkCompatiblity(self):
    oldKeys_     = ['INTERPOLATION', 'MHD_FORMULATION','RESISTIVE_MHD']
    replaceKeys_ = ['RECONSTRUCTION','DIVB_CONTROL',   'RESISTIVITY']
    if (os.path.exists(self.work_dir+'/definitions.hpp')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.hpp')

      for i in range(len(oldKeys_)):
        pf.ReplaceWord(oldKeys_[i], replaceKeys_[i])


  def GenerateOptionsList(self):
    """Creates a default option list.

    This method of class DefineProblem will create a
    default valued options list for each entry in the
    entries list. These are essentially the options
    that will be browsed in the Pluto Setup Menu.
    """
    phylist = ['HD','RHD','MHD','RMHD','ResRMHD','GRMHD']
    dimlist = ['1','2','3']
    geolist = ['CARTESIAN','CYLINDRICAL','SPHERICAL']
    mg_list = ['NO','YES']
    bfolist = ['NO','VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)']
    coolist = ['NO','POWER_LAW','TABULATED','SNEq','SYNCH'] #,'MINEq','H2_COOL', 'KROME']
    reclist = ['FLAT','LINEAR','PARABOLIC','LimO3','WENO3','WENOZ','MP5']
    if (self.flag_dict['WITH-HO']):
      tmslist = ['EULER','RK2','RK3','RK4','ARK4','IERK45']
    else:
      tmslist = ['EULER','RK2','RK3']
      
    ntrlist = ['%d'%n for n in range(9)]
    prtlist = ['NO','PARTICLES_LP','PARTICLES_CR']
    udplist = ['%d'%n for n in range(32)]
    udclist = ['%d'%n for n in range(32)]

    self.options = [phylist, dimlist, geolist, mg_list, 
                    bfolist, coolist, reclist, 
                    tmslist, ntrlist, prtlist, udplist, udclist]

  def AfterFlagLists(self):
    """Modify options and default list based on command-line flags.

    This method is called after generation of default option list.
    It modifies the members of the options list and if required
    also the default list based on the conditions required by
    the flags set using system arguments.
    """
    if self.flag_dict['FULL']:
      self.options[self.entries.index('RECONSTRUCTION')] = ['FLAT',
                                                            'LINEAR',
                                                            'PARABOLIC',
                                                            'WENO3',
                                                            'WENOZ',
                                                            'MP5']


    if self.flag_dict['WITH-FARGO']:
      self.options[self.entries.index('PHYSICS')] = ['HD', 'MHD']
      self.options[self.entries.index('DIMENSIONS')] = ['2','3']
      self.default[self.entries.index('DIMENSIONS')] = '2'

    if self.flag_dict['WITH-SB']:
      self.options[self.entries.index('PHYSICS')] = ['HD', 'MHD']
      self.default[self.entries.index('PHYSICS')] = 'MHD'
      self.options[self.entries.index('DIMENSIONS')] = ['2', '3']
      self.default[self.entries.index('DIMENSIONS')] = '2'
      self.options[self.entries.index('BODY_FORCE')] = ['VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)']
      self.default[self.entries.index('BODY_FORCE')] = 'VECTOR'


  def ReadOrBrowse(self, Ents=None, Opts=None, Defs=None, MenuTitle=None):
    """Reads or Browses the entries, options and defaults to create header file.

    This method either reads the already exsisting definitions.hpp file or browses
    throught the three lists which are provided as inputs.

    **Inputs**:
      1. Ents = List of entries. [None]
      2. Opts = List of options corresponding to each member in Ents [None]
      3. Defs = List of default value from Opts corresponding to each member in Ents [None]
      4. MenuTitle = Title for the Menu to be Browsed [None]
    """
    if (os.path.exists(self.work_dir+'/definitions.hpp')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.hpp')
      pf.UpdateListFromFile(Ents, Defs)
      for i in range(len(Ents)):
        if Defs[i] not in Opts[i]:
          Defs[i] = Opts[i][0]
        else:
          pass

    # Provides Browsing options using the menu file in case of no automatic update flag.
    if self.auto_update == 0:
      selection = ''
      menu.SetTitle(MenuTitle)
      selection = menu.Browse(Ents, default=Defs, options=Opts)

  def ProcessHDModule(self):
    """
    Provides entries, options and defaults specific to Hydro Module.
    Also updates them accordingly if required by flags.
    """
    self.entries_HD = ['DUST_FLUID','EOS', 'ENTROPY_SWITCH',
                       'THERMAL_CONDUCTION', 'VISCOSITY',
                       'RADIATION', 'ROTATING_FRAME']
    self.default_HD = ['NO','IDEAL', 'NO',
                       'NO',    'NO',
                       'NO', 'NO']
    self.options_HD = [['NO','YES'],['IDEAL','PVTE_LAW','ISOTHERMAL'],
                       ['NO','SELECTIVE','ALWAYS'],
                       ['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_LEGENDRE'],
                       ['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_LEGENDRE'],
                       ['NO','YES'],
                       ['NO','YES']]

  def ProcessRHDModule(self):
    """
    Provides entries, options and defaults specific to Relativistic
    Hydro Module. Also updates them accordingly if required by flags.
    """
    self.entries_RHD = ['EOS', 'ENTROPY_SWITCH','RADIATION']
    self.default_RHD = ['IDEAL', 'NO', 'NO']
    self.options_RHD = [['IDEAL','TAUB'],
                        ['NO','SELECTIVE','ALWAYS'],
                        ['NO','YES']]

  def ProcessMHDModule(self):
    """
    Provides entries, options and defaults specific to Magneto-
    Hydro Module.Also updates them accordingly if required by flags.
    """
    self.entries_MHD = ['EOS', 'ENTROPY_SWITCH', 'DIVB_CONTROL',
                        'BACKGROUND_FIELD',
                        'AMBIPOLAR_DIFFUSION', 'RESISTIVITY', 'HALL_MHD',
                        'THERMAL_CONDUCTION', 'VISCOSITY',
                        'RADIATION', 'ROTATING_FRAME']
    self.default_MHD = ['IDEAL','NO','EIGHT_WAVES',
                        'NO',
                        'NO','NO','NO',
                        'NO','NO',
                        'NO','NO']
    self.options_MHD = [['IDEAL','PVTE_LAW','ISOTHERMAL'],
#                            ['NO','YES'],
                        ['NO','SELECTIVE','ALWAYS'],
                        ['NO','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT'],
                        ['NO','YES'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','EXPLICIT'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','YES'],
                        ['NO','YES']]

    if self.flag_dict['WITH-SB'] or self.flag_dict['WITH-FARGO']:
      indx_ = self.entries_MHD.index('DIVB_CONTROL')
      self.options_MHD[indx_] = ['CONSTRAINED_TRANSPORT', 'NO']
      self.default_MHD[indx_] = 'CONSTRAINED_TRANSPORT'

  def ProcessRMHDModule(self):
    """
    Provides entries, options and defaults specific to Relativisitc
    Magneto-Hydro Module.Also updates them accordingly if required by flags.
    """
    self.entries_RMHD = ['EOS', 'ENTROPY_SWITCH','RADIATION','DIVB_CONTROL']
    self.default_RMHD = ['IDEAL',           'NO',       'NO',          'NO']
    self.options_RMHD = [['IDEAL', 'TAUB'],
                         ['NO','SELECTIVE','ALWAYS'],
                         ['NO','YES'],
                         ['NO','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT']]

  def ProcessResRMHDModule(self):
    """
    Provides entries, options and defaults specific to Resistive Relativistic
    Magneto-Hydro Module. Also updates them accordingly if required by flags.
    """
    self.entries_ResRMHD = ['EOS', 'ENTROPY_SWITCH','DIVB_CONTROL','DIVE_CONTROL']
    self.default_ResRMHD = ['IDEAL',           'NO',          'NO',          'NO']
    self.options_ResRMHD = [['IDEAL', 'TAUB'],
#                            ['NO','YES'],
                         ['NO','SELECTIVE','ALWAYS'],
                         ['NO','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT'],
                         ['NO','DIV_CLEANING','CONSTRAINED_TRANSPORT']]

  def ProcessGRMHDModule(self):
    """
    Provides entries, options and defaults specific to General Relativistic
    Magneto-Hydro Module.Also updates them accordingly if required by flags.
    """
    self.entries_GRMHD = ['EOS', 'ENTROPY_SWITCH','RADIATION','DIVB_CONTROL']
    self.default_GRMHD = ['IDEAL',           'NO',       'NO',          'NO']
    self.options_GRMHD = [['IDEAL', 'TAUB'],
                         ['NO','SELECTIVE','ALWAYS'],
                         ['NO','YES'],
                         ['NO','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT']]

  def ProcessParticlesModule(self):
    """
    Provides entries, options and defaults specific to the
    Particle module. Also updates them accordingly if required by flags.
    """
    if (self.default[self.entries.index('PARTICLES')] == 'PARTICLES_LP'):

      if (self.default[self.entries.index('PHYSICS')] == 'MHD' or self.default[self.entries.index('PHYSICS')] == 'RMHD'):
        self.entries_Particles = ['PARTICLES_LP_SPECTRA','PARTICLES_LP_SHOCK_SPECTRA','PARTICLES_LP_MREC_SPECTRA']
        self.default_Particles = ['YES','YES','NO']
        self.options_Particles = [['NO', 'YES'],['NO', 'YES'],['NO', 'YES']]
      else :
        self.entries_Particles = ['PARTICLES_LP_SPECTRA']
        self.default_Particles = ['NO']
        self.options_Particles = [['NO']]

    else :
      if (self.default[self.entries.index('PARTICLES')] == 'PARTICLES_CR'):
        self.entries_Particles = ['PARTICLES_CR_FEEDBACK']
        self.default_Particles = ['NO']
        self.options_Particles = [['NO', 'YES']]

  def ProcessKROMEModule(self):
    """
    Provides entries, options, and defaults for the KROME cooling module.
    """
    netwrkfiles = os.listdir(self.krome_dir+'networks/')
    self.entries_KROME = ['NETWORK_FILE', 'USE_N','COOLING_TYPE', 'HEATING_TYPE', 'GAMMA_TYPE']
    self.default_KROME = [netwrkfiles[0], 'NO', 'NONE', 'NONE', 'DEFAULT']
    self.options_KROME = [netwrkfiles, ['NO', 'YES'], [ 'NONE', 'H2', 'ATOMIC', 'Z', 'H2,ATOMIC'],
                         ['NONE', 'COMPRESS', 'PHOTO', 'XRAY', 'CHEM'],
                         ['DEFAULT','FULL', 'ROT', 'VIB', 'EXACT', 'REDUCED']]

    if self.eos == 'IDEAL':
      self.options_KROME[4] = ['DEFAULT']



  def ProcessUserDefPara(self):
    """
    Sets the Userdefined parameters
    """
    self.udef_params = ['USER_PAR_%.2d'%i for i in range(self.nparam)]
    if (os.path.exists(self.work_dir+'/definitions.hpp')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.hpp')
      scrh = pf.LocateString('parameters')
      k0   = scrh[0][0] + 2
      par_lines = pf.ReadLines(k0, k0 + self.nparam)
      for n in range(self.nparam):
        try:
          x = par_lines[n].split()
          x[0] == '#define'
        except IndexError:
          pass
        else:
          if (x[0] == "#define"):
            self.udef_params[n] = x[1]
          else:
            break;

    if self.auto_update == 0:
      menu.SetTitle ("User-defined Parameters")
      par_entries = ['%d'%i for i in range(self.nparam)]
      menu.Insert(par_entries,self.udef_params)


  def ProcessUserDefConst(self):
    """
    Sets the Userdefined Constants.
    """
    if (os.path.exists(self.work_dir+'/definitions.hpp')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.hpp')
      old_beg_scrh = pf.LocateString('symbolic')
      try:
        old_beg_scrh[0][0]
      except IndexError:
        pass
      else:
        del_indx = pf.LocateString('USER_DEF_CONSTANTS')
        pf.DeleteLines(del_indx[0][0], del_indx[0][0])
        old_beg_scrh = pf.LocateString('symbolic')
        pf.ReplaceLine('// [Beg] user-defined constants (do not change this line)', old_beg_scrh[0][0])
        old_end_scrh = pf.LocateString('supplementary')
        pf.InsertLine('// [End] user-defined constants (do not change this line)', old_end_scrh[0][0] - 1)

      scrh_beg = pf.LocateString('[Beg]')
      k_beg   = scrh_beg[0][0]+1
      scrh_end = pf.LocateString('[End]')
      k_end   = scrh_end[0][0]-1
      const_lines = pf.ReadLines(k_beg, k_end)
      #print const_lines
      for n in range(len(const_lines)):
        x = const_lines[n].split()
        try:
          x[0] == '#define'
        except IndexError:
          pass
        else:
          if (x[0] == '#define'):
            self.udef_const.append(x[1])
            self.udef_const_vals.append(x[2])
          else:
            continue


  def NonUserFriendlyConst(self):
    """
        Sets the non-user friendly constants. : OBSOLETE

        This function only stores the Non-User Friendly constants
        whose values are different than their DEFAULT values and
        incorporates them in User Defined Constants.
    """
    self.non_usfr = []
    tmplist1 = ['INITIAL_SMOOTHING', 'WARNING_MESSAGES',
                'INTERNAL_BOUNDARY', 'SHOCK_FLATTENING']
    tmplist2     = ['NO','YES','NO','NO']
    tmplist2_def = ['NO','YES','NO','NO']

    if not self.flag_dict['WITH-FD']:
      tmplist1 = tmplist1 + ['CHAR_LIMITING', 'LIMITER']
      tmplist2 = tmplist2 + ['NO', 'DEFAULT']
      tmplist2_def = tmplist2_def + ['NO', 'DEFAULT']

    if 'DIVB_CONTROL' in self.mod_entries:
      divb_mode = self.mod_default[self.mod_entries.index('DIVB_CONTROL')]
      if divb_mode == 'CONSTRAINED_TRANSPORT':
        tmplist1 = tmplist1 + ['CT_EMF_AVERAGE', 'CT_EN_CORRECTION', 'ASSIGN_VECTOR_POTENTIAL']
        tmplist2 = tmplist2 + ['UCT_HLL', 'NO', 'NO']
        tmplist2_def = tmplist2_def + ['UCT_HLL', 'NO', 'NO']
      else:
        tmplist1 = tmplist1 + ['ASSIGN_VECTOR_POTENTIAL']
        tmplist2 = tmplist2 + ['NO']
        tmplist2_def = tmplist2_def + ['NO']

      if not self.flag_dict['WITH-CHOMBO']:
        tmplist1 = tmplist1 + ['UPDATE_VECTOR_POTENTIAL']
        tmplist2 = tmplist2 + ['NO']
        tmplist2_def = tmplist2_def + ['NO']

#        if 'HANCOCK' in self.default:
#            if (self.phymodule == 'RMHD'):
#                tmplist1 = tmplist1 + ['PRIMITIVE_HANCOCK']
#                tmplist2 = tmplist2 + ['NO']
#           else:
#                tmplist1 = tmplist1 + ['PRIMITIVE_HANCOCK']
#                tmplist2 = tmplist2 + ['YES']

    longword = max(len(w) for w in tmplist1)

    if (os.path.exists(self.work_dir+'/definitions.hpp')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.hpp')
      scrh = pf.LocateString('supplementary')
      try:
        scrh[0]
      except IndexError:
        pass
      else:
        pf.UpdateListFromFile(tmplist1, tmplist2)

    for i in range(len(tmplist1)):
      if tmplist2[i] != tmplist2_def[i]:
        self.non_usfr = self.non_usfr + ['#define  '+tmplist1[i].ljust(28)+'   '+tmplist2[i]+'\n']

  def AppendAdditionalFiles(self):
    """
    Adds additional object files based on
    modular defintions and requirements.
    """
    reconstruction = self.default[self.entries.index('RECONSTRUCTION')]

    if reconstruction == 'FLAT':
      self.additional_files.append('flat_states.o')
    if reconstruction == 'LINEAR':
      self.additional_files.append('plm_states.o')
    elif reconstruction == 'PARABOLIC':
      self.additional_files.append('ppm_states.o')
    elif reconstruction in ['LimO3','WENO3','WENOZ','MP5']:
      self.additional_files.append('mp5_states.o')
    else:
      pass

    self.additional_files.append('rk_step.o')
    self.additional_files.append('update_stage.o')

    # -- Diffusion time stepping --
    parabolic_update = 0
    if 'SUPER_TIME_STEPPING' in self.mod_default:
      parabolic_update = 1
      self.additional_files += ['sts.o']

    if 'RK_LEGENDRE' in self.mod_default:
      parabolic_update = 1
      self.additional_files += ['rkl.o']

    if 'EXPLICIT' in self.mod_default:
      parabolic_update = 1

    if (parabolic_update):
      self.additional_files.append('parabolic_flux.o')


  def AppendPlutoPathAndFlags(self):
    """
    Adds additional C flags and path to 'makefile' based on
    module defintions and requirements.
    """
    self.pluto_path.append(self.phymodule+'/')

    # Add path for diffusion modules
    dis_eff = ['Dust_Fluid','Thermal_Conduction', 'Viscosity']
    for de in dis_eff:
      if (    (de.upper() in self.mod_entries)
          and (self.mod_default[self.mod_entries.index(de.upper())] != 'NO')):
        self.pluto_path.append(de+'/')

    if    self.phymodule == 'MHD' \
       or self.phymodule == 'RMHD' \
       or self.phymodule == "ResRMHD" \
       or self.phymodule == "GRMHD":
      divb_mode   = self.mod_default[self.mod_entries.index('DIVB_CONTROL')]
      if (self.phymodule == 'MHD'): \
        resistivity = self.mod_default[self.mod_entries.index('RESISTIVITY')]

      if divb_mode == 'CONSTRAINED_TRANSPORT':
        self.pluto_path.append('MHD/CT/')
      elif divb_mode == 'DIV_CLEANING':
        self.pluto_path.append('MHD/GLM/')
      else:
          pass

      if     self.phymodule == 'MHD' \
         and self.mod_default[self.mod_entries.index('AMBIPOLAR_DIFFUSION')] != 'NO':
        self.pluto_path.append('MHD/Ambipolar_Diffusion/')

      if self.phymodule == 'MHD' and \
        self.mod_default[self.mod_entries.index('HALL_MHD')] != 'NO':
        self.pluto_path.append('MHD/Hall_MHD/')

      if     self.phymodule == 'MHD' \
         and self.mod_default[self.mod_entries.index('RESISTIVITY')] != 'NO':
        self.pluto_path.append('MHD/Resistivity/')

      if self.phymodule == 'ResRMHD':
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Hot fix to reverse path list for debugging purposes.
        # This will have to be removed at some point later on...
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.pluto_path=list(reversed(self.pluto_path))

        if ('rk_step.o' in self.additional_files):
          self.additional_files.remove('rk_step.o')
          self.additional_files.append('rk_step_imex.o')


    # -- Add path for Radiation --
    if self.phymodule != "ResRMHD":
      radiation = self.mod_default[self.mod_entries.index('RADIATION')]
      if (radiation == "YES"):
        self.pluto_path.append('Radiation/')

    # -- Add path for FARGO --

    if self.flag_dict['WITH-FARGO']:
      self.additional_flags.append(' -DFARGO')
      self.pluto_path.append('Fargo/')

    # -- Add path for High Order (HO) --

    if self.flag_dict['WITH-HO']:
      self.additional_flags.append(' -DHIGH_ORDER')
      self.pluto_path.append('High_Order/')

    # -- Add path for Poisson Solver --

    if self.flag_dict['WITH-POISSON_SOLVER']:
      self.additional_flags.append(' -DPOISSON_SOLVER')
      self.pluto_path.append('PoissonSolver/')

    # -- Add path for ShearingBox (SB) --

    if self.flag_dict['WITH-SB']:
      self.pluto_path.append('MHD/ShearingBox/')
      self.additional_flags.append(' -DSHEARINGBOX')

    # -- Add path for Body Force --

    bd_force = self.default[self.entries.index('BODY_FORCE')]

    # -- Add path for Cooling --
    cool_mode = self.default[self.entries.index('COOLING')]
    if cool_mode != 'NO':
      if cool_mode not in [ "POWER_LAW", "SYNCH"]:
        self.pluto_path.append('Cooling/')
         
      if cool_mode == 'TABULATED':
        self.pluto_path.append('Cooling/Tabulated/')
      elif cool_mode == 'POWER_LAW':
        self.pluto_path.append('Cooling/Power_Law/')
      elif cool_mode == 'SYNCH':
        self.pluto_path.append('Cooling/Synch/')
      else:
        self.pluto_path.append('Cooling/'+ cool_mode +'/')

    # -- Add path for EOS --
    # if 'EOS' in self.mod_entries:
    #   if self.eos == "IDEAL":
    #     self.pluto_path.append('EOS/Ideal/')
    #   elif self.eos == "ISOTHERMAL":
    #     self.pluto_path.append('EOS/Isothermal/')
    #   elif self.eos == "TAUB":
    #     self.pluto_path.append('EOS/Taub/')

    # -- Add path for Particles --
    particles = self.default[self.entries.index('PARTICLES')]
    if (   particles == "PARTICLES_LP" or particles == "PARTICLES_CR"):
      self.pluto_path.append('Particles/')

  def UpdatePlutoIni(self):
    """
    Updates pluto.ini file with values of UserDefined Parameters
    """
    pf   = pfIO.PlutoFiles(self.work_dir+'/pluto.ini')
    scrh = pf.LocateString('[Parameters]')
    try:
      scrh[0]
    except IndexError:
      print("Parameters keyword not found in pluto.ini")
      sys.exit()
    else:
      pass

    ipos = scrh[0][0] + 2
    tmplist1 = pf.ReadLines(ipos,ipos+512)
    paradict = {}
    cmms = []
    for x in tmplist1:
      if (len(x.split()) == 0): continue  # skip blank lines
      paradict.update({x.split()[0]:x.split()[1]})
      try:
        cmmval = x.split()[2]
      except IndexError:
        cmms.append('')
        continue
      else:
        if cmmval == '#' or cmmval.startswith('#'):
          cmms.append(' '.join(x.split()[2:]))
        else:
          cmms.append('')

    for x in self.udef_params:
      if x in paradict.keys():
        pf.InsertLine(x.ljust(28) + paradict[x] +'  '+cmms[self.udef_params.index(x)]+'\n', ipos)
      else:
        try:
          cmms[self.udef_params.index(x)]
        except IndexError:
          pf.InsertLine(x.ljust(28) + '0.0' + '\n', ipos)
        else:
          pf.InsertLine(x.ljust(28) + '0.0'+'  '+cmms[self.udef_params.index(x)]+ '\n', ipos)
      ipos = ipos + 1
    pf.DeleteLines(ipos,ipos+100)


  def WriteKromeOpts(self):
    self.kromeoptstr = 'python '+self.krome_dir+'krome -interfaceC -n='+self.krome_dir+'networks/'+ self.default_KROME[self.entries_KROME.index('NETWORK_FILE')]
    if self.default_KROME[self.entries_KROME.index('USE_N')] == 'YES':
      self.kromeoptstr += ' -useN'
    else:
      self.kromeoptstr += ' -useX'

    if self.default_KROME[self.entries_KROME.index('COOLING_TYPE')] != 'NONE':
      self.kromeoptstr += ' -cooling='+ self.default_KROME[self.entries_KROME.index('COOLING_TYPE')]

    if self.default_KROME[self.entries_KROME.index('HEATING_TYPE')] != 'NONE':
      self.kromeoptstr += ' -heating='+ self.default_KROME[self.entries_KROME.index('HEATING_TYPE')]

    if self.default_KROME[self.entries_KROME.index('GAMMA_TYPE')] != 'DEFAULT':
      self.kromeoptstr += ' -gamma='+ self.default_KROME[self.entries_KROME.index('GAMMA_TYPE')]

    #remove all checks :
    self.kromeoptstr += ' -noCheck -noSinkCheck'

  def WriteDefFileList(self):
    """
    Writes all modular entries, options, defaults into a list.
    """
    for x in self.entries:
        self.def_file_list.append( '#define  '+x.ljust(28)+'   '
                                  +self.default[self.entries.index(x)]+'\n')

    self.def_file_list.append('\n// Physics dependent declarations\n\n')
    self.phymodule = self.default[self.entries.index('PHYSICS')]

    tmp1 = ['entries_%s'%self.phymodule, 'default_%s'%self.phymodule]

    # print self.phymodule
    # sys.exit(0)

    self.mod_entries = self.__getattribute__(tmp1[0])
    self.mod_default = self.__getattribute__(tmp1[1])

    for x in self.mod_entries:
      self.def_file_list.append( '#define  '+x.ljust(28)+'   '
                                + self.mod_default[self.mod_entries.index(x)]+'\n')

    self.AppendAdditionalFiles()
    self.AppendPlutoPathAndFlags()

    try:
     self.entries_KROME
    except AttributeError:
     pass
    else:
     self.def_file_list.append('\n/* -- KROME Cooling Options -- */\n\n')
     for x in self.entries_KROME:
        self.def_file_list.append('#define  '+x.ljust(28)+'   '+self.default_KROME[self.entries_KROME.index(x)]+'\n')

    try:
     self.entries_Particles
    except AttributeError:
     pass
    else:
     self.def_file_list.append('\n/* -- Particles Options -- */\n\n')
     for x in self.entries_Particles:
        self.def_file_list.append('#define  '+x.ljust(28)+'   '+self.default_Particles[self.entries_Particles.index(x)]+'\n')
 
    # always insert user-defined paramters and constants sections

    self.def_file_list.append('\n// User-defined parameters (labels)\n\n')
    for x in self.udef_params:
      self.def_file_list.append('#define  '+x.ljust(28)+'   '+'%d'%self.udef_params.index(x)+'\n')

    self.UpdatePlutoIni()


    self.def_file_list.append('\n// [Beg] user-defined constants (do not change this line)\n\n')
    self.NonUserFriendlyConst()
    if len(self.non_usfr) > 0:
      for x in self.non_usfr:
        self.def_file_list.append(x)

    for i in range(len(self.udef_const)):
      self.def_file_list.append('#define  '+self.udef_const[i].ljust(28)+'   '+self.udef_const_vals[i]+'\n')

    self.def_file_list.append('\n// [End] user-defined constants (do not change this line)\n')

    #self.def_file_list.append('\n/* -- supplementary constants (user editable) -- */ \n\n')
