from __future__ import print_function, absolute_import, division
import _pysimple
import f90wrap.runtime
import logging
import numpy

class Simple(f90wrap.runtime.FortranModule):
    """
    Module simple
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
        lines 5-163
    
    """
    @staticmethod
    def init_field(self, vmec_file, ans_s, ans_tp, amultharm, aintegmode):
        """
        init_field(self, vmec_file, ans_s, ans_tp, amultharm, aintegmode)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 26-60
        
        Parameters
        ----------
        self : Tracer
        vmec_file : str
        ans_s : int
        ans_tp : int
        amultharm : int
        aintegmode : int
        
        """
        _pysimple.f90wrap_init_field(self=self._handle, vmec_file=vmec_file, \
            ans_s=ans_s, ans_tp=ans_tp, amultharm=amultharm, aintegmode=aintegmode)
    
    @staticmethod
    def init_params(self, z_charge=None, m_mass=None, e_kin=None, npoints=None, \
        store_step=None, relerr=None):
        """
        init_params(self[, z_charge, m_mass, e_kin, npoints, store_step, relerr])
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 62-96
        
        Parameters
        ----------
        self : Tracer
        z_charge : int
        m_mass : int
        e_kin : unknown
        npoints : int
        store_step : int
        relerr : unknown
        
        """
        _pysimple.f90wrap_init_params(self=self._handle, z_charge=z_charge, \
            m_mass=m_mass, e_kin=e_kin, npoints=npoints, store_step=store_step, \
            relerr=relerr)
    
    @staticmethod
    def init_sympl(self, f, z0, dtau, dtaumin, rtol_init, mode_init):
        """
        init_sympl(self, f, z0, dtau, dtaumin, rtol_init, mode_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 98-121
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        z0 : unknown array
        dtau : unknown
        dtaumin : unknown
        rtol_init : unknown
        mode_init : int
        
        """
        _pysimple.f90wrap_init_sympl(si=self._handle, f=f._handle, z0=z0, dtau=dtau, \
            dtaumin=dtaumin, rtol_init=rtol_init, mode_init=mode_init)
    
    @staticmethod
    def init_integrator(self, z0):
        """
        init_integrator(self, z0)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 123-127
        
        Parameters
        ----------
        self : Tracer
        z0 : unknown array
        
        """
        _pysimple.f90wrap_init_integrator(self=self._handle, z0=z0)
    
    @staticmethod
    def _timestep(self, s, th, ph, lam):
        """
        ierr = _timestep(self, s, th, ph, lam)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 129-143
        
        Parameters
        ----------
        self : Tracer
        s : unknown
        th : unknown
        ph : unknown
        lam : unknown
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_timestep(self=self._handle, s=s, th=th, ph=ph, lam=lam)
        return ierr
    
    @staticmethod
    def _timestep_z(self, z):
        """
        ierr = _timestep_z(self, z)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 145-149
        
        Parameters
        ----------
        self : Tracer
        z : unknown array
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_timestep_z(self=self._handle, z=z)
        return ierr
    
    @staticmethod
    def _timestep_sympl_z(self, f, z):
        """
        ierr = _timestep_sympl_z(self, f, z)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 151-162
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        z : unknown array
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_timestep_sympl_z(si=self._handle, f=f._handle, z=z)
        return ierr
    
    @staticmethod
    def tstep(*args, **kwargs):
        """
        tstep(*args, **kwargs)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 20-23
        
        Overloaded interface containing the following procedures:
          _timestep
          _timestep_z
          _timestep_sympl_z
        
        """
        for proc in [Simple._timestep, Simple._timestep_z, Simple._timestep_sympl_z]:
            try:
                return proc(*args, **kwargs)
            except TypeError:
                continue
        
    
    _dt_array_initialisers = []
    

simple = Simple()

class Cut_Detector(f90wrap.runtime.FortranModule):
    """
    Module cut_detector
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
        lines 165-331
    
    """
    @f90wrap.runtime.register_class("pysimple.CutDetector")
    class CutDetector(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=cutdetector)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 177-183
        
        """
        def __init__(self, handle=None):
            """
            self = Cutdetector()
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                lines 177-183
            
            
            Returns
            -------
            this : Cutdetector
            	Object to be constructed
            
            
            Automatically generated constructor for cutdetector
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _pysimple.f90wrap_cutdetector_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Cutdetector
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                lines 177-183
            
            Parameters
            ----------
            this : Cutdetector
            	Object to be destructed
            
            
            Automatically generated destructor for cutdetector
            """
            if self._alloc:
                _pysimple.f90wrap_cutdetector_finalise(this=self._handle)
        
        @property
        def fper(self):
            """
            Element fper ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 178
            
            """
            return _pysimple.f90wrap_cutdetector__get__fper(self._handle)
        
        @fper.setter
        def fper(self, fper):
            _pysimple.f90wrap_cutdetector__set__fper(self._handle, fper)
        
        @property
        def alam_prev(self):
            """
            Element alam_prev ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 180
            
            """
            return _pysimple.f90wrap_cutdetector__get__alam_prev(self._handle)
        
        @alam_prev.setter
        def alam_prev(self, alam_prev):
            _pysimple.f90wrap_cutdetector__set__alam_prev(self._handle, alam_prev)
        
        @property
        def par_inv(self):
            """
            Element par_inv ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 180
            
            """
            return _pysimple.f90wrap_cutdetector__get__par_inv(self._handle)
        
        @par_inv.setter
        def par_inv(self, par_inv):
            _pysimple.f90wrap_cutdetector__set__par_inv(self._handle, par_inv)
        
        @property
        def iper(self):
            """
            Element iper ftype=integer           pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 181
            
            """
            return _pysimple.f90wrap_cutdetector__get__iper(self._handle)
        
        @iper.setter
        def iper(self, iper):
            _pysimple.f90wrap_cutdetector__set__iper(self._handle, iper)
        
        @property
        def itip(self):
            """
            Element itip ftype=integer           pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 181
            
            """
            return _pysimple.f90wrap_cutdetector__get__itip(self._handle)
        
        @itip.setter
        def itip(self, itip):
            _pysimple.f90wrap_cutdetector__set__itip(self._handle, itip)
        
        @property
        def kper(self):
            """
            Element kper ftype=integer           pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 181
            
            """
            return _pysimple.f90wrap_cutdetector__get__kper(self._handle)
        
        @kper.setter
        def kper(self, kper):
            _pysimple.f90wrap_cutdetector__set__kper(self._handle, kper)
        
        @property
        def orb_sten(self):
            """
            Element orb_sten ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 182
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_cutdetector__array__orb_sten(self._handle)
            if array_handle in self._arrays:
                orb_sten = self._arrays[array_handle]
            else:
                orb_sten = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_cutdetector__array__orb_sten)
                self._arrays[array_handle] = orb_sten
            return orb_sten
        
        @orb_sten.setter
        def orb_sten(self, orb_sten):
            self.orb_sten[...] = orb_sten
        
        @property
        def coef(self):
            """
            Element coef ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 182
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_cutdetector__array__coef(self._handle)
            if array_handle in self._arrays:
                coef = self._arrays[array_handle]
            else:
                coef = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_cutdetector__array__coef)
                self._arrays[array_handle] = coef
            return coef
        
        @coef.setter
        def coef(self, coef):
            self.coef[...] = coef
        
        @property
        def ipoi(self):
            """
            Element ipoi ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
                line 183
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_cutdetector__array__ipoi(self._handle)
            if array_handle in self._arrays:
                ipoi = self._arrays[array_handle]
            else:
                ipoi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_cutdetector__array__ipoi)
                self._arrays[array_handle] = ipoi
            return ipoi
        
        @ipoi.setter
        def ipoi(self, ipoi):
            self.ipoi[...] = ipoi
        
        def __str__(self):
            ret = ['<cutdetector>{\n']
            ret.append('    fper : ')
            ret.append(repr(self.fper))
            ret.append(',\n    alam_prev : ')
            ret.append(repr(self.alam_prev))
            ret.append(',\n    par_inv : ')
            ret.append(repr(self.par_inv))
            ret.append(',\n    iper : ')
            ret.append(repr(self.iper))
            ret.append(',\n    itip : ')
            ret.append(repr(self.itip))
            ret.append(',\n    kper : ')
            ret.append(repr(self.kper))
            ret.append(',\n    orb_sten : ')
            ret.append(repr(self.orb_sten))
            ret.append(',\n    coef : ')
            ret.append(repr(self.coef))
            ret.append(',\n    ipoi : ')
            ret.append(repr(self.ipoi))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def init(self, fper, z):
        """
        init(self, fper, z)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 186-211
        
        Parameters
        ----------
        self : Cutdetector
        fper : unknown
        z : unknown array
        
        ---------------------------------------------------------------------------
         Prepare calculation of orbit tip by interpolation and buffer for Poincare plot:
        """
        _pysimple.f90wrap_init(self=self._handle, fper=fper, z=z)
    
    @staticmethod
    def trace_to_cut(self, si, f, z, var_cut):
        """
        cut_type, ierr = trace_to_cut(self, si, f, z, var_cut)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 213-277
        
        Parameters
        ----------
        self : Cutdetector
        si : Symplecticintegrator
        f : Fieldcan
        z : unknown array
        var_cut : unknown array
        
        Returns
        -------
        cut_type : int
        ierr : int
        
        -------------------------------------------------------------------------
         Tip detection and interpolation
        """
        cut_type, ierr = _pysimple.f90wrap_trace_to_cut(self=self._handle, \
            si=si._handle, f=f._handle, z=z, var_cut=var_cut)
        return cut_type, ierr
    
    @staticmethod
    def fract_dimension(ntr, rt, fraction):
        """
        fract_dimension(ntr, rt, fraction)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 279-330
        
        Parameters
        ----------
        ntr : int
        rt : unknown array
        fraction : unknown
        
        """
        _pysimple.f90wrap_fract_dimension(ntr=ntr, rt=rt, fraction=fraction)
    
    @property
    def n_tip_vars(self):
        """
        Element n_tip_vars ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            line 173
        
        """
        return _pysimple.f90wrap_cut_detector__get__n_tip_vars()
    
    @property
    def nplagr(self):
        """
        Element nplagr ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            line 174
        
        """
        return _pysimple.f90wrap_cut_detector__get__nplagr()
    
    @property
    def nder(self):
        """
        Element nder ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            line 175
        
        """
        return _pysimple.f90wrap_cut_detector__get__nder()
    
    def __str__(self):
        ret = ['<cut_detector>{\n']
        ret.append('    n_tip_vars : ')
        ret.append(repr(self.n_tip_vars))
        ret.append(',\n    nplagr : ')
        ret.append(repr(self.nplagr))
        ret.append(',\n    nder : ')
        ret.append(repr(self.nder))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

cut_detector = Cut_Detector()

class Simple_Main(f90wrap.runtime.FortranModule):
    """
    Module simple_main
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
        lines 333-1001
    
    """
    @staticmethod
    def run(self):
        """
        run(self)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 348-391
        
        Parameters
        ----------
        norb : Tracer
        
        """
        _pysimple.f90wrap_run(norb=self._handle)
    
    @staticmethod
    def finalize():
        """
        finalize()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 393-398
        
        
        """
        _pysimple.f90wrap_finalize()
    
    @staticmethod
    def write_output():
        """
        write_output()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 400-415
        
        
        """
        _pysimple.f90wrap_write_output()
    
    @staticmethod
    def init_starting_surf():
        """
        init_starting_surf()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 417-431
        
        
        """
        _pysimple.f90wrap_init_starting_surf()
    
    @staticmethod
    def init_starting_points_ants(unit):
        """
        init_starting_points_ants(unit)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 433-453
        
        Parameters
        ----------
        unit : int
        
        """
        _pysimple.f90wrap_init_starting_points_ants(unit=unit)
    
    @staticmethod
    def init_starting_points():
        """
        init_starting_points()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 455-508
        
        
        """
        _pysimple.f90wrap_init_starting_points()
    
    @staticmethod
    def init_starting_points_global():
        """
        init_starting_points_global()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 510-576
        
        
        """
        _pysimple.f90wrap_init_starting_points_global()
    
    @staticmethod
    def trace_orbit(self, ipart):
        """
        trace_orbit(self, ipart)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
            lines 578-1000
        
        Parameters
        ----------
        anorb : Tracer
        ipart : int
        
        --------------------------------
         Initialize tip detector
        """
        _pysimple.f90wrap_trace_orbit(anorb=self._handle, ipart=ipart)
    
    _dt_array_initialisers = []
    

simple_main = Simple_Main()

class Orbit_Symplectic(f90wrap.runtime.FortranModule):
    """
    Module orbit_symplectic
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
        lines 5-1481
    
    """
    @f90wrap.runtime.register_class("pysimple.SymplecticIntegrator")
    class SymplecticIntegrator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=symplecticintegrator)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 13-34
        
        """
        def __init__(self, handle=None):
            """
            self = Symplecticintegrator()
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                lines 13-34
            
            
            Returns
            -------
            this : Symplecticintegrator
            	Object to be constructed
            
            
            Automatically generated constructor for symplecticintegrator
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _pysimple.f90wrap_symplecticintegrator_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Symplecticintegrator
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                lines 13-34
            
            Parameters
            ----------
            this : Symplecticintegrator
            	Object to be destructed
            
            
            Automatically generated destructor for symplecticintegrator
            """
            if self._alloc:
                _pysimple.f90wrap_symplecticintegrator_finalise(this=self._handle)
        
        @property
        def nlag(self):
            """
            Element nlag ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 14
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__nlag(self._handle)
        
        @nlag.setter
        def nlag(self, nlag):
            _pysimple.f90wrap_symplecticintegrator__set__nlag(self._handle, nlag)
        
        @property
        def nbuf(self):
            """
            Element nbuf ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 15
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__nbuf(self._handle)
        
        @nbuf.setter
        def nbuf(self, nbuf):
            _pysimple.f90wrap_symplecticintegrator__set__nbuf(self._handle, nbuf)
        
        @property
        def extrap_field(self):
            """
            Element extrap_field ftype=logical pytype=bool
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 16
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__extrap_field(self._handle)
        
        @extrap_field.setter
        def extrap_field(self, extrap_field):
            _pysimple.f90wrap_symplecticintegrator__set__extrap_field(self._handle, \
                extrap_field)
        
        @property
        def atol(self):
            """
            Element atol ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 17
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__atol(self._handle)
        
        @atol.setter
        def atol(self, atol):
            _pysimple.f90wrap_symplecticintegrator__set__atol(self._handle, atol)
        
        @property
        def rtol(self):
            """
            Element rtol ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 18
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__rtol(self._handle)
        
        @rtol.setter
        def rtol(self, rtol):
            _pysimple.f90wrap_symplecticintegrator__set__rtol(self._handle, rtol)
        
        @property
        def z(self):
            """
            Element z ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_symplecticintegrator__array__z(self._handle)
            if array_handle in self._arrays:
                z = self._arrays[array_handle]
            else:
                z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_symplecticintegrator__array__z)
                self._arrays[array_handle] = z
            return z
        
        @z.setter
        def z(self, z):
            self.z[...] = z
        
        @property
        def pthold(self):
            """
            Element pthold ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 21
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__pthold(self._handle)
        
        @pthold.setter
        def pthold(self, pthold):
            _pysimple.f90wrap_symplecticintegrator__set__pthold(self._handle, pthold)
        
        @property
        def kbuf(self):
            """
            Element kbuf ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 23
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__kbuf(self._handle)
        
        @kbuf.setter
        def kbuf(self, kbuf):
            _pysimple.f90wrap_symplecticintegrator__set__kbuf(self._handle, kbuf)
        
        @property
        def kt(self):
            """
            Element kt ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 24
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__kt(self._handle)
        
        @kt.setter
        def kt(self, kt):
            _pysimple.f90wrap_symplecticintegrator__set__kt(self._handle, kt)
        
        @property
        def k(self):
            """
            Element k ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 25
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__k(self._handle)
        
        @k.setter
        def k(self, k):
            _pysimple.f90wrap_symplecticintegrator__set__k(self._handle, k)
        
        @property
        def bufind(self):
            """
            Element bufind ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_symplecticintegrator__array__bufind(self._handle)
            if array_handle in self._arrays:
                bufind = self._arrays[array_handle]
            else:
                bufind = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_symplecticintegrator__array__bufind)
                self._arrays[array_handle] = bufind
            return bufind
        
        @bufind.setter
        def bufind(self, bufind):
            self.bufind[...] = bufind
        
        @property
        def zbuf(self):
            """
            Element zbuf ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 27
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_symplecticintegrator__array__zbuf(self._handle)
            if array_handle in self._arrays:
                zbuf = self._arrays[array_handle]
            else:
                zbuf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_symplecticintegrator__array__zbuf)
                self._arrays[array_handle] = zbuf
            return zbuf
        
        @zbuf.setter
        def zbuf(self, zbuf):
            self.zbuf[...] = zbuf
        
        @property
        def coef(self):
            """
            Element coef ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_symplecticintegrator__array__coef(self._handle)
            if array_handle in self._arrays:
                coef = self._arrays[array_handle]
            else:
                coef = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_symplecticintegrator__array__coef)
                self._arrays[array_handle] = coef
            return coef
        
        @coef.setter
        def coef(self, coef):
            self.coef[...] = coef
        
        @property
        def ntau(self):
            """
            Element ntau ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 30
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__ntau(self._handle)
        
        @ntau.setter
        def ntau(self, ntau):
            _pysimple.f90wrap_symplecticintegrator__set__ntau(self._handle, ntau)
        
        @property
        def dt(self):
            """
            Element dt ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 31
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__dt(self._handle)
        
        @dt.setter
        def dt(self, dt):
            _pysimple.f90wrap_symplecticintegrator__set__dt(self._handle, dt)
        
        @property
        def pabs(self):
            """
            Element pabs ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 32
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__pabs(self._handle)
        
        @pabs.setter
        def pabs(self, pabs):
            _pysimple.f90wrap_symplecticintegrator__set__pabs(self._handle, pabs)
        
        @property
        def mode(self):
            """
            Element mode ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 34
            
            """
            return _pysimple.f90wrap_symplecticintegrator__get__mode(self._handle)
        
        @mode.setter
        def mode(self, mode):
            _pysimple.f90wrap_symplecticintegrator__set__mode(self._handle, mode)
        
        def __str__(self):
            ret = ['<symplecticintegrator>{\n']
            ret.append('    nlag : ')
            ret.append(repr(self.nlag))
            ret.append(',\n    nbuf : ')
            ret.append(repr(self.nbuf))
            ret.append(',\n    extrap_field : ')
            ret.append(repr(self.extrap_field))
            ret.append(',\n    atol : ')
            ret.append(repr(self.atol))
            ret.append(',\n    rtol : ')
            ret.append(repr(self.rtol))
            ret.append(',\n    z : ')
            ret.append(repr(self.z))
            ret.append(',\n    pthold : ')
            ret.append(repr(self.pthold))
            ret.append(',\n    kbuf : ')
            ret.append(repr(self.kbuf))
            ret.append(',\n    kt : ')
            ret.append(repr(self.kt))
            ret.append(',\n    k : ')
            ret.append(repr(self.k))
            ret.append(',\n    bufind : ')
            ret.append(repr(self.bufind))
            ret.append(',\n    zbuf : ')
            ret.append(repr(self.zbuf))
            ret.append(',\n    coef : ')
            ret.append(repr(self.coef))
            ret.append(',\n    ntau : ')
            ret.append(repr(self.ntau))
            ret.append(',\n    dt : ')
            ret.append(repr(self.dt))
            ret.append(',\n    pabs : ')
            ret.append(repr(self.pabs))
            ret.append(',\n    mode : ')
            ret.append(repr(self.mode))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("pysimple.MultistageIntegrator")
    class MultistageIntegrator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=multistageintegrator)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 41-44
        
        """
        def __init__(self, handle=None):
            """
            self = Multistageintegrator()
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                lines 41-44
            
            
            Returns
            -------
            this : Multistageintegrator
            	Object to be constructed
            
            
            Automatically generated constructor for multistageintegrator
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _pysimple.f90wrap_multistageintegrator_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Multistageintegrator
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                lines 41-44
            
            Parameters
            ----------
            this : Multistageintegrator
            	Object to be destructed
            
            
            Automatically generated destructor for multistageintegrator
            """
            if self._alloc:
                _pysimple.f90wrap_multistageintegrator_finalise(this=self._handle)
        
        @property
        def s(self):
            """
            Element s ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 42
            
            """
            return _pysimple.f90wrap_multistageintegrator__get__s(self._handle)
        
        @s.setter
        def s(self, s):
            _pysimple.f90wrap_multistageintegrator__set__s(self._handle, s)
        
        @property
        def alpha(self):
            """
            Element alpha ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_multistageintegrator__array__alpha(self._handle)
            if array_handle in self._arrays:
                alpha = self._arrays[array_handle]
            else:
                alpha = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_multistageintegrator__array__alpha)
                self._arrays[array_handle] = alpha
            return alpha
        
        @alpha.setter
        def alpha(self, alpha):
            self.alpha[...] = alpha
        
        @property
        def beta(self):
            """
            Element beta ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 43
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_multistageintegrator__array__beta(self._handle)
            if array_handle in self._arrays:
                beta = self._arrays[array_handle]
            else:
                beta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_multistageintegrator__array__beta)
                self._arrays[array_handle] = beta
            return beta
        
        @beta.setter
        def beta(self, beta):
            self.beta[...] = beta
        
        def init_array_stages(self):
            self.stages = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _pysimple.f90wrap_multistageintegrator__array_getitem__stages,
                                            _pysimple.f90wrap_multistageintegrator__array_setitem__stages,
                                            _pysimple.f90wrap_multistageintegrator__array_len__stages,
                                            """
            Element stages ftype=type(symplecticintegrator) pytype=Symplecticintegrator
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
                line 44
            
            """, Orbit_Symplectic.SymplecticIntegrator)
            return self.stages
        
        def __str__(self):
            ret = ['<multistageintegrator>{\n']
            ret.append('    s : ')
            ret.append(repr(self.s))
            ret.append(',\n    alpha : ')
            ret.append(repr(self.alpha))
            ret.append(',\n    beta : ')
            ret.append(repr(self.beta))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_stages]
        
    
    @staticmethod
    def orbit_sympl_init(self, f, z, dt, ntau, rtol_init, mode_init, nlag):
        """
        orbit_sympl_init(self, f, z, dt, ntau, rtol_init, mode_init, nlag)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 49-76
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        z : unknown array
        dt : unknown
        ntau : int
        rtol_init : unknown
        mode_init : int
        nlag : int
        
        """
        _pysimple.f90wrap_orbit_sympl_init(si=self._handle, f=f._handle, z=z, dt=dt, \
            ntau=ntau, rtol_init=rtol_init, mode_init=mode_init, nlag=nlag)
    
    @staticmethod
    def f_sympl_euler1(self, f, n, x, fvec, iflag):
        """
        f_sympl_euler1(self, f, n, x, fvec, iflag)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 80-91
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        n : int
        x : unknown array
        fvec : unknown array
        iflag : int
        
        """
        _pysimple.f90wrap_f_sympl_euler1(si=self._handle, f=f._handle, n=n, x=x, \
            fvec=fvec, iflag=iflag)
    
    @staticmethod
    def jac_sympl_euler1(self, f, x, jac):
        """
        jac_sympl_euler1(self, f, x, jac)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 95-108
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        x : unknown array
        jac : unknown array
        
        """
        _pysimple.f90wrap_jac_sympl_euler1(si=self._handle, f=f._handle, x=x, jac=jac)
    
    @staticmethod
    def f_sympl_euler2(self, f, n, x, fvec, iflag):
        """
        f_sympl_euler2(self, f, n, x, fvec, iflag)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 112-124
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        n : int
        x : unknown array
        fvec : unknown array
        iflag : int
        
        """
        _pysimple.f90wrap_f_sympl_euler2(si=self._handle, f=f._handle, n=n, x=x, \
            fvec=fvec, iflag=iflag)
    
    @staticmethod
    def jac_sympl_euler2(self, f, x, jac):
        """
        jac_sympl_euler2(self, f, x, jac)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 128-139
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        x : unknown array
        jac : unknown array
        
        """
        _pysimple.f90wrap_jac_sympl_euler2(si=self._handle, f=f._handle, x=x, jac=jac)
    
    @staticmethod
    def f_midpoint_part1(self, f, n, x, fvec):
        """
        f_midpoint_part1(self, f, n, x, fvec)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 143-156
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        n : int
        x : unknown array
        fvec : unknown array
        
        """
        _pysimple.f90wrap_f_midpoint_part1(si=self._handle, f=f._handle, n=n, x=x, \
            fvec=fvec)
    
    @staticmethod
    def f_midpoint_part2(self, f, n, x, fvec):
        """
        f_midpoint_part2(self, f, n, x, fvec)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 160-174
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        n : int
        x : unknown array
        fvec : unknown array
        
        """
        _pysimple.f90wrap_f_midpoint_part2(si=self._handle, f=f._handle, n=n, x=x, \
            fvec=fvec)
    
    @staticmethod
    def jac_midpoint_part1(self, f, x, jac):
        """
        jac_midpoint_part1(self, f, x, jac)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 178-225
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        x : unknown array
        jac : unknown array
        
        """
        _pysimple.f90wrap_jac_midpoint_part1(si=self._handle, f=f._handle, x=x, jac=jac)
    
    @staticmethod
    def jac_midpoint_part2(self, f, fmid, x, jac):
        """
        jac_midpoint_part2(self, f, fmid, x, jac)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 229-249
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        fmid : Fieldcan
        x : unknown array
        jac : unknown array
        
        """
        _pysimple.f90wrap_jac_midpoint_part2(si=self._handle, f=f._handle, \
            fmid=fmid._handle, x=x, jac=jac)
    
    @staticmethod
    def newton1(self, f, x, maxit, xlast):
        """
        newton1(self, f, x, maxit, xlast)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 253-296
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        x : unknown array
        maxit : int
        xlast : unknown array
        
        """
        _pysimple.f90wrap_newton1(si=self._handle, f=f._handle, x=x, maxit=maxit, \
            xlast=xlast)
    
    @staticmethod
    def newton2(self, f, x, atol, rtol, maxit, xlast):
        """
        newton2(self, f, x, atol, rtol, maxit, xlast)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 298-341
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        x : unknown array
        atol : unknown
        rtol : unknown
        maxit : int
        xlast : unknown array
        
        """
        _pysimple.f90wrap_newton2(si=self._handle, f=f._handle, x=x, atol=atol, \
            rtol=rtol, maxit=maxit, xlast=xlast)
    
    @staticmethod
    def newton_midpoint(self, f, x, atol, rtol, maxit, xlast):
        """
        newton_midpoint(self, f, x, atol, rtol, maxit, xlast)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 343-378
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        x : unknown array
        atol : unknown
        rtol : unknown
        maxit : int
        xlast : unknown array
        
        """
        _pysimple.f90wrap_newton_midpoint(si=self._handle, f=f._handle, x=x, atol=atol, \
            rtol=rtol, maxit=maxit, xlast=xlast)
    
    @staticmethod
    def coeff_rk_gauss(n, a, b, c):
        """
        coeff_rk_gauss(n, a, b, c)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 380-442
        
        Parameters
        ----------
        n : int
        a : unknown array
        b : unknown array
        c : unknown array
        
        """
        _pysimple.f90wrap_coeff_rk_gauss(n=n, a=a, b=b, c=c)
    
    @staticmethod
    def coeff_rk_lobatto(n, a, ahat, b, c):
        """
        coeff_rk_lobatto(n, a, ahat, b, c)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 645-678
        
        Parameters
        ----------
        n : int
        a : unknown array
        ahat : unknown array
        b : unknown array
        c : unknown array
        
        """
        _pysimple.f90wrap_coeff_rk_lobatto(n=n, a=a, ahat=ahat, b=b, c=c)
    
    @staticmethod
    def orbit_timestep_sympl(self, f):
        """
        ierr = orbit_timestep_sympl(self, f)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 889-914
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl(si=self._handle, f=f._handle)
        return ierr
    
    @staticmethod
    def orbit_timestep_sympl_multi(self, f):
        """
        ierr = orbit_timestep_sympl_multi(self, f)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 918-931
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl_multi(mi=self._handle, \
            f=f._handle)
        return ierr
    
    @staticmethod
    def orbit_sympl_init_multi(self, f, z, dtau, ntau, rtol_init, alpha, beta):
        """
        orbit_sympl_init_multi(self, f, z, dtau, ntau, rtol_init, alpha, beta)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 935-951
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        alpha : unknown array
        beta : unknown array
        
        """
        _pysimple.f90wrap_orbit_sympl_init_multi(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init, alpha=alpha, beta=beta)
    
    @staticmethod
    def orbit_sympl_init_verlet(self, f, z, dtau, ntau, rtol_init):
        """
        orbit_sympl_init_verlet(self, f, z, dtau, ntau, rtol_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 955-966
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        
        """
        _pysimple.f90wrap_orbit_sympl_init_verlet(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init)
    
    @staticmethod
    def orbit_sympl_init_order4(self, f, z, dtau, ntau, rtol_init):
        """
        orbit_sympl_init_order4(self, f, z, dtau, ntau, rtol_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 970-986
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        
        """
        _pysimple.f90wrap_orbit_sympl_init_order4(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init)
    
    @staticmethod
    def orbit_sympl_init_mclachlan4(self, f, z, dtau, ntau, rtol_init):
        """
        orbit_sympl_init_mclachlan4(self, f, z, dtau, ntau, rtol_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 990-1012
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        
        """
        _pysimple.f90wrap_orbit_sympl_init_mclachlan4(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init)
    
    @staticmethod
    def orbit_sympl_init_blanes4(self, f, z, dtau, ntau, rtol_init):
        """
        orbit_sympl_init_blanes4(self, f, z, dtau, ntau, rtol_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1016-1041
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        
        """
        _pysimple.f90wrap_orbit_sympl_init_blanes4(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init)
    
    @staticmethod
    def orbit_sympl_init_kahan6(self, f, z, dtau, ntau, rtol_init):
        """
        orbit_sympl_init_kahan6(self, f, z, dtau, ntau, rtol_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1045-1067
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        
        """
        _pysimple.f90wrap_orbit_sympl_init_kahan6(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init)
    
    @staticmethod
    def orbit_sympl_init_kahan8(self, f, z, dtau, ntau, rtol_init):
        """
        orbit_sympl_init_kahan8(self, f, z, dtau, ntau, rtol_init)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1071-1101
        
        Parameters
        ----------
        mi : Multistageintegrator
        f : Fieldcan
        z : unknown array
        dtau : unknown
        ntau : int
        rtol_init : unknown
        
        """
        _pysimple.f90wrap_orbit_sympl_init_kahan8(mi=self._handle, f=f._handle, z=z, \
            dtau=dtau, ntau=ntau, rtol_init=rtol_init)
    
    @staticmethod
    def orbit_timestep_sympl_euler1(self, f):
        """
        ierr = orbit_timestep_sympl_euler1(self, f)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1105-1167
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl_euler1(si=self._handle, \
            f=f._handle)
        return ierr
    
    @staticmethod
    def orbit_timestep_sympl_euler2(self, f):
        """
        ierr = orbit_timestep_sympl_euler2(self, f)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1171-1230
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl_euler2(si=self._handle, \
            f=f._handle)
        return ierr
    
    @staticmethod
    def orbit_timestep_sympl_midpoint(self, f):
        """
        ierr = orbit_timestep_sympl_midpoint(self, f)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1234-1290
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl_midpoint(si=self._handle, \
            f=f._handle)
        return ierr
    
    @staticmethod
    def orbit_timestep_sympl_rk_gauss(self, f, s):
        """
        ierr = orbit_timestep_sympl_rk_gauss(self, f, s)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1294-1417
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        s : int
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl_rk_gauss(si=self._handle, \
            f=f._handle, s=s)
        return ierr
    
    @staticmethod
    def orbit_timestep_sympl_rk_lobatto(self, f, s):
        """
        ierr = orbit_timestep_sympl_rk_lobatto(self, f, s)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            lines 1421-1465
        
        Parameters
        ----------
        si : Symplecticintegrator
        f : Fieldcan
        s : int
        
        Returns
        -------
        ierr : int
        
        """
        ierr = _pysimple.f90wrap_orbit_timestep_sympl_rk_lobatto(si=self._handle, \
            f=f._handle, s=s)
        return ierr
    
    @property
    def nlag_max(self):
        """
        Element nlag_max ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            line 11
        
        """
        return _pysimple.f90wrap_orbit_symplectic__get__nlag_max()
    
    @property
    def nbuf_max(self):
        """
        Element nbuf_max ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            line 12
        
        """
        return _pysimple.f90wrap_orbit_symplectic__get__nbuf_max()
    
    @property
    def s_max(self):
        """
        Element s_max ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/orbit_symplectic.f90.i \
            line 40
        
        """
        return _pysimple.f90wrap_orbit_symplectic__get__s_max()
    
    def __str__(self):
        ret = ['<orbit_symplectic>{\n']
        ret.append('    nlag_max : ')
        ret.append(repr(self.nlag_max))
        ret.append(',\n    nbuf_max : ')
        ret.append(repr(self.nbuf_max))
        ret.append(',\n    s_max : ')
        ret.append(repr(self.s_max))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

orbit_symplectic = Orbit_Symplectic()

class Field_Can_Mod(f90wrap.runtime.FortranModule):
    """
    Module field_can_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
        lines 5-371
    
    """
    @f90wrap.runtime.register_class("pysimple.FieldCan")
    class FieldCan(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=fieldcan)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 8-26
        
        """
        def __init__(self, handle=None):
            """
            self = Fieldcan()
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                lines 8-26
            
            
            Returns
            -------
            this : Fieldcan
            	Object to be constructed
            
            
            Automatically generated constructor for fieldcan
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _pysimple.f90wrap_fieldcan_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Fieldcan
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                lines 8-26
            
            Parameters
            ----------
            this : Fieldcan
            	Object to be destructed
            
            
            Automatically generated destructor for fieldcan
            """
            if self._alloc:
                _pysimple.f90wrap_fieldcan_finalise(this=self._handle)
        
        @property
        def field_type(self):
            """
            Element field_type ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 9
            
            """
            return _pysimple.f90wrap_fieldcan__get__field_type(self._handle)
        
        @field_type.setter
        def field_type(self, field_type):
            _pysimple.f90wrap_fieldcan__set__field_type(self._handle, field_type)
        
        @property
        def ath(self):
            """
            Element ath ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 10
            
            """
            return _pysimple.f90wrap_fieldcan__get__ath(self._handle)
        
        @ath.setter
        def ath(self, ath):
            _pysimple.f90wrap_fieldcan__set__ath(self._handle, ath)
        
        @property
        def aph(self):
            """
            Element aph ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 10
            
            """
            return _pysimple.f90wrap_fieldcan__get__aph(self._handle)
        
        @aph.setter
        def aph(self, aph):
            _pysimple.f90wrap_fieldcan__set__aph(self._handle, aph)
        
        @property
        def hth(self):
            """
            Element hth ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 11
            
            """
            return _pysimple.f90wrap_fieldcan__get__hth(self._handle)
        
        @hth.setter
        def hth(self, hth):
            _pysimple.f90wrap_fieldcan__set__hth(self._handle, hth)
        
        @property
        def hph(self):
            """
            Element hph ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 11
            
            """
            return _pysimple.f90wrap_fieldcan__get__hph(self._handle)
        
        @hph.setter
        def hph(self, hph):
            _pysimple.f90wrap_fieldcan__set__hph(self._handle, hph)
        
        @property
        def bmod(self):
            """
            Element bmod ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 12
            
            """
            return _pysimple.f90wrap_fieldcan__get__bmod(self._handle)
        
        @bmod.setter
        def bmod(self, bmod):
            _pysimple.f90wrap_fieldcan__set__bmod(self._handle, bmod)
        
        @property
        def dath(self):
            """
            Element dath ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dath(self._handle)
            if array_handle in self._arrays:
                dath = self._arrays[array_handle]
            else:
                dath = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dath)
                self._arrays[array_handle] = dath
            return dath
        
        @dath.setter
        def dath(self, dath):
            self.dath[...] = dath
        
        @property
        def daph(self):
            """
            Element daph ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 13
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__daph(self._handle)
            if array_handle in self._arrays:
                daph = self._arrays[array_handle]
            else:
                daph = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__daph)
                self._arrays[array_handle] = daph
            return daph
        
        @daph.setter
        def daph(self, daph):
            self.daph[...] = daph
        
        @property
        def dhth(self):
            """
            Element dhth ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dhth(self._handle)
            if array_handle in self._arrays:
                dhth = self._arrays[array_handle]
            else:
                dhth = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dhth)
                self._arrays[array_handle] = dhth
            return dhth
        
        @dhth.setter
        def dhth(self, dhth):
            self.dhth[...] = dhth
        
        @property
        def dhph(self):
            """
            Element dhph ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dhph(self._handle)
            if array_handle in self._arrays:
                dhph = self._arrays[array_handle]
            else:
                dhph = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dhph)
                self._arrays[array_handle] = dhph
            return dhph
        
        @dhph.setter
        def dhph(self, dhph):
            self.dhph[...] = dhph
        
        @property
        def dbmod(self):
            """
            Element dbmod ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dbmod(self._handle)
            if array_handle in self._arrays:
                dbmod = self._arrays[array_handle]
            else:
                dbmod = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dbmod)
                self._arrays[array_handle] = dbmod
            return dbmod
        
        @dbmod.setter
        def dbmod(self, dbmod):
            self.dbmod[...] = dbmod
        
        @property
        def d2ath(self):
            """
            Element d2ath ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2ath(self._handle)
            if array_handle in self._arrays:
                d2ath = self._arrays[array_handle]
            else:
                d2ath = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2ath)
                self._arrays[array_handle] = d2ath
            return d2ath
        
        @d2ath.setter
        def d2ath(self, d2ath):
            self.d2ath[...] = d2ath
        
        @property
        def d2aph(self):
            """
            Element d2aph ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2aph(self._handle)
            if array_handle in self._arrays:
                d2aph = self._arrays[array_handle]
            else:
                d2aph = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2aph)
                self._arrays[array_handle] = d2aph
            return d2aph
        
        @d2aph.setter
        def d2aph(self, d2aph):
            self.d2aph[...] = d2aph
        
        @property
        def d2hth(self):
            """
            Element d2hth ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2hth(self._handle)
            if array_handle in self._arrays:
                d2hth = self._arrays[array_handle]
            else:
                d2hth = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2hth)
                self._arrays[array_handle] = d2hth
            return d2hth
        
        @d2hth.setter
        def d2hth(self, d2hth):
            self.d2hth[...] = d2hth
        
        @property
        def d2hph(self):
            """
            Element d2hph ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2hph(self._handle)
            if array_handle in self._arrays:
                d2hph = self._arrays[array_handle]
            else:
                d2hph = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2hph)
                self._arrays[array_handle] = d2hph
            return d2hph
        
        @d2hph.setter
        def d2hph(self, d2hph):
            self.d2hph[...] = d2hph
        
        @property
        def d2bmod(self):
            """
            Element d2bmod ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2bmod(self._handle)
            if array_handle in self._arrays:
                d2bmod = self._arrays[array_handle]
            else:
                d2bmod = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2bmod)
                self._arrays[array_handle] = d2bmod
            return d2bmod
        
        @d2bmod.setter
        def d2bmod(self, d2bmod):
            self.d2bmod[...] = d2bmod
        
        @property
        def h(self):
            """
            Element h ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 20
            
            """
            return _pysimple.f90wrap_fieldcan__get__h(self._handle)
        
        @h.setter
        def h(self, h):
            _pysimple.f90wrap_fieldcan__set__h(self._handle, h)
        
        @property
        def pth(self):
            """
            Element pth ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 20
            
            """
            return _pysimple.f90wrap_fieldcan__get__pth(self._handle)
        
        @pth.setter
        def pth(self, pth):
            _pysimple.f90wrap_fieldcan__set__pth(self._handle, pth)
        
        @property
        def vpar(self):
            """
            Element vpar ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 20
            
            """
            return _pysimple.f90wrap_fieldcan__get__vpar(self._handle)
        
        @vpar.setter
        def vpar(self, vpar):
            _pysimple.f90wrap_fieldcan__set__vpar(self._handle, vpar)
        
        @property
        def dvpar(self):
            """
            Element dvpar ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dvpar(self._handle)
            if array_handle in self._arrays:
                dvpar = self._arrays[array_handle]
            else:
                dvpar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dvpar)
                self._arrays[array_handle] = dvpar
            return dvpar
        
        @dvpar.setter
        def dvpar(self, dvpar):
            self.dvpar[...] = dvpar
        
        @property
        def dh(self):
            """
            Element dh ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dh(self._handle)
            if array_handle in self._arrays:
                dh = self._arrays[array_handle]
            else:
                dh = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dh)
                self._arrays[array_handle] = dh
            return dh
        
        @dh.setter
        def dh(self, dh):
            self.dh[...] = dh
        
        @property
        def dpth(self):
            """
            Element dpth ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__dpth(self._handle)
            if array_handle in self._arrays:
                dpth = self._arrays[array_handle]
            else:
                dpth = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__dpth)
                self._arrays[array_handle] = dpth
            return dpth
        
        @dpth.setter
        def dpth(self, dpth):
            self.dpth[...] = dpth
        
        @property
        def d2vpar(self):
            """
            Element d2vpar ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2vpar(self._handle)
            if array_handle in self._arrays:
                d2vpar = self._arrays[array_handle]
            else:
                d2vpar = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2vpar)
                self._arrays[array_handle] = d2vpar
            return d2vpar
        
        @d2vpar.setter
        def d2vpar(self, d2vpar):
            self.d2vpar[...] = d2vpar
        
        @property
        def d2h(self):
            """
            Element d2h ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2h(self._handle)
            if array_handle in self._arrays:
                d2h = self._arrays[array_handle]
            else:
                d2h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2h)
                self._arrays[array_handle] = d2h
            return d2h
        
        @d2h.setter
        def d2h(self, d2h):
            self.d2h[...] = d2h
        
        @property
        def d2pth(self):
            """
            Element d2pth ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _pysimple.f90wrap_fieldcan__array__d2pth(self._handle)
            if array_handle in self._arrays:
                d2pth = self._arrays[array_handle]
            else:
                d2pth = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _pysimple.f90wrap_fieldcan__array__d2pth)
                self._arrays[array_handle] = d2pth
            return d2pth
        
        @d2pth.setter
        def d2pth(self, d2pth):
            self.d2pth[...] = d2pth
        
        @property
        def mu(self):
            """
            Element mu ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 26
            
            """
            return _pysimple.f90wrap_fieldcan__get__mu(self._handle)
        
        @mu.setter
        def mu(self, mu):
            _pysimple.f90wrap_fieldcan__set__mu(self._handle, mu)
        
        @property
        def ro0(self):
            """
            Element ro0 ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
                line 26
            
            """
            return _pysimple.f90wrap_fieldcan__get__ro0(self._handle)
        
        @ro0.setter
        def ro0(self, ro0):
            _pysimple.f90wrap_fieldcan__set__ro0(self._handle, ro0)
        
        def __str__(self):
            ret = ['<fieldcan>{\n']
            ret.append('    field_type : ')
            ret.append(repr(self.field_type))
            ret.append(',\n    ath : ')
            ret.append(repr(self.ath))
            ret.append(',\n    aph : ')
            ret.append(repr(self.aph))
            ret.append(',\n    hth : ')
            ret.append(repr(self.hth))
            ret.append(',\n    hph : ')
            ret.append(repr(self.hph))
            ret.append(',\n    bmod : ')
            ret.append(repr(self.bmod))
            ret.append(',\n    dath : ')
            ret.append(repr(self.dath))
            ret.append(',\n    daph : ')
            ret.append(repr(self.daph))
            ret.append(',\n    dhth : ')
            ret.append(repr(self.dhth))
            ret.append(',\n    dhph : ')
            ret.append(repr(self.dhph))
            ret.append(',\n    dbmod : ')
            ret.append(repr(self.dbmod))
            ret.append(',\n    d2ath : ')
            ret.append(repr(self.d2ath))
            ret.append(',\n    d2aph : ')
            ret.append(repr(self.d2aph))
            ret.append(',\n    d2hth : ')
            ret.append(repr(self.d2hth))
            ret.append(',\n    d2hph : ')
            ret.append(repr(self.d2hph))
            ret.append(',\n    d2bmod : ')
            ret.append(repr(self.d2bmod))
            ret.append(',\n    h : ')
            ret.append(repr(self.h))
            ret.append(',\n    pth : ')
            ret.append(repr(self.pth))
            ret.append(',\n    vpar : ')
            ret.append(repr(self.vpar))
            ret.append(',\n    dvpar : ')
            ret.append(repr(self.dvpar))
            ret.append(',\n    dh : ')
            ret.append(repr(self.dh))
            ret.append(',\n    dpth : ')
            ret.append(repr(self.dpth))
            ret.append(',\n    d2vpar : ')
            ret.append(repr(self.d2vpar))
            ret.append(',\n    d2h : ')
            ret.append(repr(self.d2h))
            ret.append(',\n    d2pth : ')
            ret.append(repr(self.d2pth))
            ret.append(',\n    mu : ')
            ret.append(repr(self.mu))
            ret.append(',\n    ro0 : ')
            ret.append(repr(self.ro0))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def fieldcan_init(self, mu=None, ro0=None, vpar=None, field_type=None):
        """
        fieldcan_init(self[, mu, ro0, vpar, field_type])
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 29-52
        
        Parameters
        ----------
        f : Fieldcan
        mu : unknown
        ro0 : unknown
        vpar : unknown
        field_type : int
        
        """
        _pysimple.f90wrap_fieldcan_init(f=self._handle, mu=mu, ro0=ro0, vpar=vpar, \
            field_type=field_type)
    
    @staticmethod
    def get_val(self, pphi):
        """
        get_val(self, pphi)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 56-65
        
        Parameters
        ----------
        f : Fieldcan
        pphi : unknown
        
        """
        _pysimple.f90wrap_get_val(f=self._handle, pphi=pphi)
    
    @staticmethod
    def get_derivatives(self, pphi):
        """
        get_derivatives(self, pphi)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 69-82
        
        Parameters
        ----------
        f : Fieldcan
        pphi : unknown
        
        """
        _pysimple.f90wrap_get_derivatives(f=self._handle, pphi=pphi)
    
    @staticmethod
    def get_derivatives2(self, pphi):
        """
        get_derivatives2(self, pphi)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 86-120
        
        Parameters
        ----------
        f : Fieldcan
        pphi : unknown
        
        """
        _pysimple.f90wrap_get_derivatives2(f=self._handle, pphi=pphi)
    
    @staticmethod
    def eval_field_can(self, r, th_c, ph_c, mode_secders):
        """
        eval_field_can(self, r, th_c, ph_c, mode_secders)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 124-214
        
        Parameters
        ----------
        f : Fieldcan
        r : unknown
        th_c : unknown
        ph_c : unknown
        mode_secders : int
        
        """
        _pysimple.f90wrap_eval_field_can(f=self._handle, r=r, th_c=th_c, ph_c=ph_c, \
            mode_secders=mode_secders)
    
    @staticmethod
    def eval_field_booz(self, r, th_c, ph_c, mode_secders):
        """
        eval_field_booz(self, r, th_c, ph_c, mode_secders)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 218-279
        
        Parameters
        ----------
        f : Fieldcan
        r : unknown
        th_c : unknown
        ph_c : unknown
        mode_secders : int
        
        """
        _pysimple.f90wrap_eval_field_booz(f=self._handle, r=r, th_c=th_c, ph_c=ph_c, \
            mode_secders=mode_secders)
    
    @staticmethod
    def eval_field_test(self, r, th, ph, mode_secders):
        """
        eval_field_test(self, r, th, ph, mode_secders)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 282-342
        
        Parameters
        ----------
        f : Fieldcan
        r : unknown
        th : unknown
        ph : unknown
        mode_secders : int
        
        """
        _pysimple.f90wrap_eval_field_test(f=self._handle, r=r, th=th, ph=ph, \
            mode_secders=mode_secders)
    
    @staticmethod
    def eval_field(self, r, th_c, ph_c, mode_secders):
        """
        eval_field(self, r, th_c, ph_c, mode_secders)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/field_can.f90.i \
            lines 346-371
        
        Parameters
        ----------
        f : Fieldcan
        r : unknown
        th_c : unknown
        ph_c : unknown
        mode_secders : int
        
        """
        _pysimple.f90wrap_eval_field(f=self._handle, r=r, th_c=th_c, ph_c=ph_c, \
            mode_secders=mode_secders)
    
    _dt_array_initialisers = []
    

field_can_mod = Field_Can_Mod()

class Exchange_Get_Cancoord_Mod(f90wrap.runtime.FortranModule):
    """
    Module exchange_get_cancoord_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 6-12
    
    """
    @property
    def onlytheta(self):
        """
        Element onlytheta ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 8
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__onlytheta()
    
    @onlytheta.setter
    def onlytheta(self, onlytheta):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__onlytheta(onlytheta)
    
    @property
    def vartheta_c(self):
        """
        Element vartheta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__vartheta_c()
    
    @vartheta_c.setter
    def vartheta_c(self, vartheta_c):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__vartheta_c(vartheta_c)
    
    @property
    def varphi_c(self):
        """
        Element varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__varphi_c()
    
    @varphi_c.setter
    def varphi_c(self, varphi_c):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__varphi_c(varphi_c)
    
    @property
    def sqg(self):
        """
        Element sqg ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__sqg()
    
    @sqg.setter
    def sqg(self, sqg):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__sqg(sqg)
    
    @property
    def aiota(self):
        """
        Element aiota ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__aiota()
    
    @aiota.setter
    def aiota(self, aiota):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__aiota(aiota)
    
    @property
    def bcovar_vartheta(self):
        """
        Element bcovar_vartheta ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__bcovar_vartheta()
    
    @bcovar_vartheta.setter
    def bcovar_vartheta(self, bcovar_vartheta):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__bcovar_vartheta(bcovar_vartheta)
    
    @property
    def bcovar_varphi(self):
        """
        Element bcovar_varphi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__bcovar_varphi()
    
    @bcovar_varphi.setter
    def bcovar_varphi(self, bcovar_varphi):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__bcovar_varphi(bcovar_varphi)
    
    @property
    def a_theta(self):
        """
        Element a_theta ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__a_theta()
    
    @a_theta.setter
    def a_theta(self, a_theta):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__a_theta(a_theta)
    
    @property
    def a_phi(self):
        """
        Element a_phi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__a_phi()
    
    @a_phi.setter
    def a_phi(self, a_phi):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__a_phi(a_phi)
    
    @property
    def theta(self):
        """
        Element theta ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__theta()
    
    @theta.setter
    def theta(self, theta):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__theta(theta)
    
    @property
    def bctrvr_vartheta(self):
        """
        Element bctrvr_vartheta ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__bctrvr_vartheta()
    
    @bctrvr_vartheta.setter
    def bctrvr_vartheta(self, bctrvr_vartheta):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__bctrvr_vartheta(bctrvr_vartheta)
    
    @property
    def bctrvr_varphi(self):
        """
        Element bctrvr_varphi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
            line 10
        
        """
        return _pysimple.f90wrap_exchange_get_cancoord_mod__get__bctrvr_varphi()
    
    @bctrvr_varphi.setter
    def bctrvr_varphi(self, bctrvr_varphi):
        _pysimple.f90wrap_exchange_get_cancoord_mod__set__bctrvr_varphi(bctrvr_varphi)
    
    def __str__(self):
        ret = ['<exchange_get_cancoord_mod>{\n']
        ret.append('    onlytheta : ')
        ret.append(repr(self.onlytheta))
        ret.append(',\n    vartheta_c : ')
        ret.append(repr(self.vartheta_c))
        ret.append(',\n    varphi_c : ')
        ret.append(repr(self.varphi_c))
        ret.append(',\n    sqg : ')
        ret.append(repr(self.sqg))
        ret.append(',\n    aiota : ')
        ret.append(repr(self.aiota))
        ret.append(',\n    bcovar_vartheta : ')
        ret.append(repr(self.bcovar_vartheta))
        ret.append(',\n    bcovar_varphi : ')
        ret.append(repr(self.bcovar_varphi))
        ret.append(',\n    a_theta : ')
        ret.append(repr(self.a_theta))
        ret.append(',\n    a_phi : ')
        ret.append(repr(self.a_phi))
        ret.append(',\n    theta : ')
        ret.append(repr(self.theta))
        ret.append(',\n    bctrvr_vartheta : ')
        ret.append(repr(self.bctrvr_vartheta))
        ret.append(',\n    bctrvr_varphi : ')
        ret.append(repr(self.bctrvr_varphi))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

exchange_get_cancoord_mod = Exchange_Get_Cancoord_Mod()

class Params(f90wrap.runtime.FortranModule):
    """
    Module params
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
        lines 5-126
    
    """
    @f90wrap.runtime.register_class("pysimple.Tracer")
    class Tracer(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=tracer)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            lines 15-23
        
        """
        def __init__(self, handle=None):
            """
            self = Tracer()
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                lines 15-23
            
            
            Returns
            -------
            this : Tracer
            	Object to be constructed
            
            
            Automatically generated constructor for tracer
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _pysimple.f90wrap_tracer_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Tracer
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                lines 15-23
            
            Parameters
            ----------
            this : Tracer
            	Object to be destructed
            
            
            Automatically generated destructor for tracer
            """
            if self._alloc:
                _pysimple.f90wrap_tracer_finalise(this=self._handle)
        
        @property
        def fper(self):
            """
            Element fper ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 16
            
            """
            return _pysimple.f90wrap_tracer__get__fper(self._handle)
        
        @fper.setter
        def fper(self, fper):
            _pysimple.f90wrap_tracer__set__fper(self._handle, fper)
        
        @property
        def dtau(self):
            """
            Element dtau ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 17
            
            """
            return _pysimple.f90wrap_tracer__get__dtau(self._handle)
        
        @dtau.setter
        def dtau(self, dtau):
            _pysimple.f90wrap_tracer__set__dtau(self._handle, dtau)
        
        @property
        def dtaumin(self):
            """
            Element dtaumin ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 17
            
            """
            return _pysimple.f90wrap_tracer__get__dtaumin(self._handle)
        
        @dtaumin.setter
        def dtaumin(self, dtaumin):
            _pysimple.f90wrap_tracer__set__dtaumin(self._handle, dtaumin)
        
        @property
        def v0(self):
            """
            Element v0 ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 17
            
            """
            return _pysimple.f90wrap_tracer__get__v0(self._handle)
        
        @v0.setter
        def v0(self, v0):
            _pysimple.f90wrap_tracer__set__v0(self._handle, v0)
        
        @property
        def n_e(self):
            """
            Element n_e ftype=integer           pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 18
            
            """
            return _pysimple.f90wrap_tracer__get__n_e(self._handle)
        
        @n_e.setter
        def n_e(self, n_e):
            _pysimple.f90wrap_tracer__set__n_e(self._handle, n_e)
        
        @property
        def n_d(self):
            """
            Element n_d ftype=integer           pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 18
            
            """
            return _pysimple.f90wrap_tracer__get__n_d(self._handle)
        
        @n_d.setter
        def n_d(self, n_d):
            _pysimple.f90wrap_tracer__set__n_d(self._handle, n_d)
        
        @property
        def integmode(self):
            """
            Element integmode ftype=integer  pytype=int
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 19
            
            """
            return _pysimple.f90wrap_tracer__get__integmode(self._handle)
        
        @integmode.setter
        def integmode(self, integmode):
            _pysimple.f90wrap_tracer__set__integmode(self._handle, integmode)
        
        @property
        def relerr(self):
            """
            Element relerr ftype=double precision pytype=unknown
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 20
            
            """
            return _pysimple.f90wrap_tracer__get__relerr(self._handle)
        
        @relerr.setter
        def relerr(self, relerr):
            _pysimple.f90wrap_tracer__set__relerr(self._handle, relerr)
        
        @property
        def f(self):
            """
            Element f ftype=type(fieldcan) pytype=Fieldcan
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 21
            
            """
            f_handle = _pysimple.f90wrap_tracer__get__f(self._handle)
            if tuple(f_handle) in self._objs:
                f = self._objs[tuple(f_handle)]
            else:
                f = field_can_mod.FieldCan.from_handle(f_handle)
                self._objs[tuple(f_handle)] = f
            return f
        
        @f.setter
        def f(self, f):
            f = f._handle
            _pysimple.f90wrap_tracer__set__f(self._handle, f)
        
        @property
        def si(self):
            """
            Element si ftype=type(symplecticintegrator) pytype=Symplecticintegrator
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 22
            
            """
            si_handle = _pysimple.f90wrap_tracer__get__si(self._handle)
            if tuple(si_handle) in self._objs:
                si = self._objs[tuple(si_handle)]
            else:
                si = orbit_symplectic.SymplecticIntegrator.from_handle(si_handle)
                self._objs[tuple(si_handle)] = si
            return si
        
        @si.setter
        def si(self, si):
            si = si._handle
            _pysimple.f90wrap_tracer__set__si(self._handle, si)
        
        @property
        def mi(self):
            """
            Element mi ftype=type(multistageintegrator) pytype=Multistageintegrator
            
            
            Defined at \
                /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
                line 23
            
            """
            mi_handle = _pysimple.f90wrap_tracer__get__mi(self._handle)
            if tuple(mi_handle) in self._objs:
                mi = self._objs[tuple(mi_handle)]
            else:
                mi = orbit_symplectic.MultistageIntegrator.from_handle(mi_handle)
                self._objs[tuple(mi_handle)] = mi
            return mi
        
        @mi.setter
        def mi(self, mi):
            mi = mi._handle
            _pysimple.f90wrap_tracer__set__mi(self._handle, mi)
        
        def __str__(self):
            ret = ['<tracer>{\n']
            ret.append('    fper : ')
            ret.append(repr(self.fper))
            ret.append(',\n    dtau : ')
            ret.append(repr(self.dtau))
            ret.append(',\n    dtaumin : ')
            ret.append(repr(self.dtaumin))
            ret.append(',\n    v0 : ')
            ret.append(repr(self.v0))
            ret.append(',\n    n_e : ')
            ret.append(repr(self.n_e))
            ret.append(',\n    n_d : ')
            ret.append(repr(self.n_d))
            ret.append(',\n    integmode : ')
            ret.append(repr(self.integmode))
            ret.append(',\n    relerr : ')
            ret.append(repr(self.relerr))
            ret.append(',\n    f : ')
            ret.append(repr(self.f))
            ret.append(',\n    si : ')
            ret.append(repr(self.si))
            ret.append(',\n    mi : ')
            ret.append(repr(self.mi))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def read_config():
        """
        read_config()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            lines 77-89
        
        
        """
        _pysimple.f90wrap_read_config()
    
    @staticmethod
    def params_init():
        """
        params_init()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            lines 91-126
        
        
        """
        _pysimple.f90wrap_params_init()
    
    @property
    def npoi(self):
        """
        Element npoi ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 25
        
        """
        return _pysimple.f90wrap_params__get__npoi()
    
    @npoi.setter
    def npoi(self, npoi):
        _pysimple.f90wrap_params__set__npoi(npoi)
    
    @property
    def l1i(self):
        """
        Element l1i ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 25
        
        """
        return _pysimple.f90wrap_params__get__l1i()
    
    @l1i.setter
    def l1i(self, l1i):
        _pysimple.f90wrap_params__set__l1i(l1i)
    
    @property
    def nper(self):
        """
        Element nper ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 25
        
        """
        return _pysimple.f90wrap_params__get__nper()
    
    @nper.setter
    def nper(self, nper):
        _pysimple.f90wrap_params__set__nper(nper)
    
    @property
    def i(self):
        """
        Element i ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 25
        
        """
        return _pysimple.f90wrap_params__get__i()
    
    @i.setter
    def i(self, i):
        _pysimple.f90wrap_params__set__i(i)
    
    @property
    def ntestpart(self):
        """
        Element ntestpart ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 25
        
        """
        return _pysimple.f90wrap_params__get__ntestpart()
    
    @ntestpart.setter
    def ntestpart(self, ntestpart):
        _pysimple.f90wrap_params__set__ntestpart(ntestpart)
    
    @property
    def loopskip(self):
        """
        Element loopskip ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 26
        
        """
        return _pysimple.f90wrap_params__get__loopskip()
    
    @loopskip.setter
    def loopskip(self, loopskip):
        _pysimple.f90wrap_params__set__loopskip(loopskip)
    
    @property
    def iskip(self):
        """
        Element iskip ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 26
        
        """
        return _pysimple.f90wrap_params__get__iskip()
    
    @iskip.setter
    def iskip(self, iskip):
        _pysimple.f90wrap_params__set__iskip(iskip)
    
    @property
    def dphi(self):
        """
        Element dphi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_params__get__dphi()
    
    @dphi.setter
    def dphi(self, dphi):
        _pysimple.f90wrap_params__set__dphi(dphi)
    
    @property
    def phibeg(self):
        """
        Element phibeg ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_params__get__phibeg()
    
    @phibeg.setter
    def phibeg(self, phibeg):
        _pysimple.f90wrap_params__set__phibeg(phibeg)
    
    @property
    def bmod00(self):
        """
        Element bmod00 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_params__get__bmod00()
    
    @bmod00.setter
    def bmod00(self, bmod00):
        _pysimple.f90wrap_params__set__bmod00(bmod00)
    
    @property
    def rlarm(self):
        """
        Element rlarm ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_params__get__rlarm()
    
    @rlarm.setter
    def rlarm(self, rlarm):
        _pysimple.f90wrap_params__set__rlarm(rlarm)
    
    @property
    def bmax(self):
        """
        Element bmax ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_params__get__bmax()
    
    @bmax.setter
    def bmax(self, bmax):
        _pysimple.f90wrap_params__set__bmax(bmax)
    
    @property
    def bmin(self):
        """
        Element bmin ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_params__get__bmin()
    
    @bmin.setter
    def bmin(self, bmin):
        _pysimple.f90wrap_params__set__bmin(bmin)
    
    @property
    def tau(self):
        """
        Element tau ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 28
        
        """
        return _pysimple.f90wrap_params__get__tau()
    
    @tau.setter
    def tau(self, tau):
        _pysimple.f90wrap_params__set__tau(tau)
    
    @property
    def dtau(self):
        """
        Element dtau ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 28
        
        """
        return _pysimple.f90wrap_params__get__dtau()
    
    @dtau.setter
    def dtau(self, dtau):
        _pysimple.f90wrap_params__set__dtau(dtau)
    
    @property
    def dtaumin(self):
        """
        Element dtaumin ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 28
        
        """
        return _pysimple.f90wrap_params__get__dtaumin()
    
    @dtaumin.setter
    def dtaumin(self, dtaumin):
        _pysimple.f90wrap_params__set__dtaumin(dtaumin)
    
    @property
    def xi(self):
        """
        Element xi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 28
        
        """
        return _pysimple.f90wrap_params__get__xi()
    
    @xi.setter
    def xi(self, xi):
        _pysimple.f90wrap_params__set__xi(xi)
    
    @property
    def rt0(self):
        """
        Element rt0 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_params__get__rt0()
    
    @rt0.setter
    def rt0(self, rt0):
        _pysimple.f90wrap_params__set__rt0(rt0)
    
    @property
    def r0i(self):
        """
        Element r0i ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_params__get__r0i()
    
    @r0i.setter
    def r0i(self, r0i):
        _pysimple.f90wrap_params__set__r0i(r0i)
    
    @property
    def cbfi(self):
        """
        Element cbfi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_params__get__cbfi()
    
    @cbfi.setter
    def cbfi(self, cbfi):
        _pysimple.f90wrap_params__set__cbfi(cbfi)
    
    @property
    def bz0i(self):
        """
        Element bz0i ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_params__get__bz0i()
    
    @bz0i.setter
    def bz0i(self, bz0i):
        _pysimple.f90wrap_params__set__bz0i(bz0i)
    
    @property
    def bf0(self):
        """
        Element bf0 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_params__get__bf0()
    
    @bf0.setter
    def bf0(self, bf0):
        _pysimple.f90wrap_params__set__bf0(bf0)
    
    @property
    def rbig(self):
        """
        Element rbig ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_params__get__rbig()
    
    @rbig.setter
    def rbig(self, rbig):
        _pysimple.f90wrap_params__set__rbig(rbig)
    
    @property
    def sbeg(self):
        """
        Element sbeg ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 30
        
        """
        return _pysimple.f90wrap_params__get__sbeg()
    
    @sbeg.setter
    def sbeg(self, sbeg):
        _pysimple.f90wrap_params__set__sbeg(sbeg)
    
    @property
    def thetabeg(self):
        """
        Element thetabeg ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 30
        
        """
        return _pysimple.f90wrap_params__get__thetabeg()
    
    @thetabeg.setter
    def thetabeg(self, thetabeg):
        _pysimple.f90wrap_params__set__thetabeg(thetabeg)
    
    @property
    def bstart(self):
        """
        Element bstart ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__bstart(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bstart = self._arrays[array_handle]
        else:
            bstart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__bstart)
            self._arrays[array_handle] = bstart
        return bstart
    
    @bstart.setter
    def bstart(self, bstart):
        self.bstart[...] = bstart
    
    @property
    def volstart(self):
        """
        Element volstart ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 31
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__volstart(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            volstart = self._arrays[array_handle]
        else:
            volstart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__volstart)
            self._arrays[array_handle] = volstart
        return volstart
    
    @volstart.setter
    def volstart(self, volstart):
        self.volstart[...] = volstart
    
    @property
    def xstart(self):
        """
        Element xstart ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 32
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__xstart(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            xstart = self._arrays[array_handle]
        else:
            xstart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__xstart)
            self._arrays[array_handle] = xstart
        return xstart
    
    @xstart.setter
    def xstart(self, xstart):
        self.xstart[...] = xstart
    
    @property
    def zstart(self):
        """
        Element zstart ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 33
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__zstart(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zstart = self._arrays[array_handle]
        else:
            zstart = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__zstart)
            self._arrays[array_handle] = zstart
        return zstart
    
    @zstart.setter
    def zstart(self, zstart):
        self.zstart[...] = zstart
    
    @property
    def zend(self):
        """
        Element zend ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 33
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__zend(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zend = self._arrays[array_handle]
        else:
            zend = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__zend)
            self._arrays[array_handle] = zend
        return zend
    
    @zend.setter
    def zend(self, zend):
        self.zend[...] = zend
    
    @property
    def confpart_trap(self):
        """
        Element confpart_trap ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__confpart_trap(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            confpart_trap = self._arrays[array_handle]
        else:
            confpart_trap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__confpart_trap)
            self._arrays[array_handle] = confpart_trap
        return confpart_trap
    
    @confpart_trap.setter
    def confpart_trap(self, confpart_trap):
        self.confpart_trap[...] = confpart_trap
    
    @property
    def confpart_pass(self):
        """
        Element confpart_pass ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__confpart_pass(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            confpart_pass = self._arrays[array_handle]
        else:
            confpart_pass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__confpart_pass)
            self._arrays[array_handle] = confpart_pass
        return confpart_pass
    
    @confpart_pass.setter
    def confpart_pass(self, confpart_pass):
        self.confpart_pass[...] = confpart_pass
    
    @property
    def times_lost(self):
        """
        Element times_lost ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__times_lost(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            times_lost = self._arrays[array_handle]
        else:
            times_lost = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__times_lost)
            self._arrays[array_handle] = times_lost
        return times_lost
    
    @times_lost.setter
    def times_lost(self, times_lost):
        self.times_lost[...] = times_lost
    
    @property
    def contr_pp(self):
        """
        Element contr_pp ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 36
        
        """
        return _pysimple.f90wrap_params__get__contr_pp()
    
    @contr_pp.setter
    def contr_pp(self, contr_pp):
        _pysimple.f90wrap_params__set__contr_pp(contr_pp)
    
    @property
    def ibins(self):
        """
        Element ibins ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 37
        
        """
        return _pysimple.f90wrap_params__get__ibins()
    
    @ibins.setter
    def ibins(self, ibins):
        _pysimple.f90wrap_params__set__ibins(ibins)
    
    @property
    def startmode(self):
        """
        Element startmode ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 38
        
        """
        return _pysimple.f90wrap_params__get__startmode()
    
    @startmode.setter
    def startmode(self, startmode):
        _pysimple.f90wrap_params__set__startmode(startmode)
    
    @property
    def ntau(self):
        """
        Element ntau ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 39
        
        """
        return _pysimple.f90wrap_params__get__ntau()
    
    @ntau.setter
    def ntau(self, ntau):
        _pysimple.f90wrap_params__set__ntau(ntau)
    
    @property
    def integmode(self):
        """
        Element integmode ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 40
        
        """
        return _pysimple.f90wrap_params__get__integmode()
    
    @integmode.setter
    def integmode(self, integmode):
        _pysimple.f90wrap_params__set__integmode(integmode)
    
    @property
    def kpart(self):
        """
        Element kpart ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 41
        
        """
        return _pysimple.f90wrap_params__get__kpart()
    
    @kpart.setter
    def kpart(self, kpart):
        _pysimple.f90wrap_params__set__kpart(kpart)
    
    @property
    def relerr(self):
        """
        Element relerr ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 42
        
        """
        return _pysimple.f90wrap_params__get__relerr()
    
    @relerr.setter
    def relerr(self, relerr):
        _pysimple.f90wrap_params__set__relerr(relerr)
    
    @property
    def trap_par(self):
        """
        Element trap_par ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__trap_par(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            trap_par = self._arrays[array_handle]
        else:
            trap_par = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__trap_par)
            self._arrays[array_handle] = trap_par
        return trap_par
    
    @trap_par.setter
    def trap_par(self, trap_par):
        self.trap_par[...] = trap_par
    
    @property
    def perp_inv(self):
        """
        Element perp_inv ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 43
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__perp_inv(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            perp_inv = self._arrays[array_handle]
        else:
            perp_inv = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__perp_inv)
            self._arrays[array_handle] = perp_inv
        return perp_inv
    
    @perp_inv.setter
    def perp_inv(self, perp_inv):
        self.perp_inv[...] = perp_inv
    
    @property
    def iclass(self):
        """
        Element iclass ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 44
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_params__array__iclass(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            iclass = self._arrays[array_handle]
        else:
            iclass = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_params__array__iclass)
            self._arrays[array_handle] = iclass
        return iclass
    
    @iclass.setter
    def iclass(self, iclass):
        self.iclass[...] = iclass
    
    @property
    def n_tip_vars(self):
        """
        Element n_tip_vars ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 45
        
        """
        return _pysimple.f90wrap_params__get__n_tip_vars()
    
    @property
    def nplagr(self):
        """
        Element nplagr ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 46
        
        """
        return _pysimple.f90wrap_params__get__nplagr()
    
    @nplagr.setter
    def nplagr(self, nplagr):
        _pysimple.f90wrap_params__set__nplagr(nplagr)
    
    @property
    def nder(self):
        """
        Element nder ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 46
        
        """
        return _pysimple.f90wrap_params__get__nder()
    
    @nder.setter
    def nder(self, nder):
        _pysimple.f90wrap_params__set__nder(nder)
    
    @property
    def npl_half(self):
        """
        Element npl_half ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 46
        
        """
        return _pysimple.f90wrap_params__get__npl_half()
    
    @npl_half.setter
    def npl_half(self, npl_half):
        _pysimple.f90wrap_params__set__npl_half(npl_half)
    
    @property
    def norbper(self):
        """
        Element norbper ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 47
        
        """
        return _pysimple.f90wrap_params__get__norbper()
    
    @norbper.setter
    def norbper(self, norbper):
        _pysimple.f90wrap_params__set__norbper(norbper)
    
    @property
    def nfp(self):
        """
        Element nfp ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 47
        
        """
        return _pysimple.f90wrap_params__get__nfp()
    
    @nfp.setter
    def nfp(self, nfp):
        _pysimple.f90wrap_params__set__nfp(nfp)
    
    @property
    def fper(self):
        """
        Element fper ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 48
        
        """
        return _pysimple.f90wrap_params__get__fper()
    
    @fper.setter
    def fper(self, fper):
        _pysimple.f90wrap_params__set__fper(fper)
    
    @property
    def zerolam(self):
        """
        Element zerolam ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 48
        
        """
        return _pysimple.f90wrap_params__get__zerolam()
    
    @zerolam.setter
    def zerolam(self, zerolam):
        _pysimple.f90wrap_params__set__zerolam(zerolam)
    
    @property
    def tcut(self):
        """
        Element tcut ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 49
        
        """
        return _pysimple.f90wrap_params__get__tcut()
    
    @tcut.setter
    def tcut(self, tcut):
        _pysimple.f90wrap_params__set__tcut(tcut)
    
    @property
    def ntcut(self):
        """
        Element ntcut ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 50
        
        """
        return _pysimple.f90wrap_params__get__ntcut()
    
    @ntcut.setter
    def ntcut(self, ntcut):
        _pysimple.f90wrap_params__set__ntcut(ntcut)
    
    @property
    def class_plot(self):
        """
        Element class_plot ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 51
        
        """
        return _pysimple.f90wrap_params__get__class_plot()
    
    @class_plot.setter
    def class_plot(self, class_plot):
        _pysimple.f90wrap_params__set__class_plot(class_plot)
    
    @property
    def cut_in_per(self):
        """
        Element cut_in_per ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 52
        
        """
        return _pysimple.f90wrap_params__get__cut_in_per()
    
    @cut_in_per.setter
    def cut_in_per(self, cut_in_per):
        _pysimple.f90wrap_params__set__cut_in_per(cut_in_per)
    
    @property
    def local(self):
        """
        Element local ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 53
        
        """
        return _pysimple.f90wrap_params__get__local()
    
    @local.setter
    def local(self, local):
        _pysimple.f90wrap_params__set__local(local)
    
    @property
    def fast_class(self):
        """
        Element fast_class ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 54
        
        """
        return _pysimple.f90wrap_params__get__fast_class()
    
    @fast_class.setter
    def fast_class(self, fast_class):
        _pysimple.f90wrap_params__set__fast_class(fast_class)
    
    @property
    def swcoll(self):
        """
        Element swcoll ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 56
        
        """
        return _pysimple.f90wrap_params__get__swcoll()
    
    @swcoll.setter
    def swcoll(self, swcoll):
        _pysimple.f90wrap_params__set__swcoll(swcoll)
    
    @property
    def am1(self):
        """
        Element am1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__am1()
    
    @property
    def am2(self):
        """
        Element am2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__am2()
    
    @property
    def z1(self):
        """
        Element z1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__z1()
    
    @property
    def z2(self):
        """
        Element z2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__z2()
    
    @property
    def densi1(self):
        """
        Element densi1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__densi1()
    
    @property
    def densi2(self):
        """
        Element densi2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__densi2()
    
    @property
    def tempi1(self):
        """
        Element tempi1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__tempi1()
    
    @property
    def tempi2(self):
        """
        Element tempi2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__tempi2()
    
    @property
    def tempe(self):
        """
        Element tempe ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_params__get__tempe()
    
    @property
    def dchichi(self):
        """
        Element dchichi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_params__get__dchichi()
    
    @dchichi.setter
    def dchichi(self, dchichi):
        _pysimple.f90wrap_params__set__dchichi(dchichi)
    
    @property
    def slowrate(self):
        """
        Element slowrate ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_params__get__slowrate()
    
    @slowrate.setter
    def slowrate(self, slowrate):
        _pysimple.f90wrap_params__set__slowrate(slowrate)
    
    @property
    def dchichi_norm(self):
        """
        Element dchichi_norm ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_params__get__dchichi_norm()
    
    @dchichi_norm.setter
    def dchichi_norm(self, dchichi_norm):
        _pysimple.f90wrap_params__set__dchichi_norm(dchichi_norm)
    
    @property
    def slowrate_norm(self):
        """
        Element slowrate_norm ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_params__get__slowrate_norm()
    
    @slowrate_norm.setter
    def slowrate_norm(self, slowrate_norm):
        _pysimple.f90wrap_params__set__slowrate_norm(slowrate_norm)
    
    @property
    def deterministic(self):
        """
        Element deterministic ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 60
        
        """
        return _pysimple.f90wrap_params__get__deterministic()
    
    @deterministic.setter
    def deterministic(self, deterministic):
        _pysimple.f90wrap_params__set__deterministic(deterministic)
    
    @property
    def notrace_passing(self):
        """
        Element notrace_passing ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 62
        
        """
        return _pysimple.f90wrap_params__get__notrace_passing()
    
    @notrace_passing.setter
    def notrace_passing(self, notrace_passing):
        _pysimple.f90wrap_params__set__notrace_passing(notrace_passing)
    
    @property
    def face_al(self):
        """
        Element face_al ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 63
        
        """
        return _pysimple.f90wrap_params__get__face_al()
    
    @face_al.setter
    def face_al(self, face_al):
        _pysimple.f90wrap_params__set__face_al(face_al)
    
    @property
    def trace_time(self):
        """
        Element trace_time ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 63
        
        """
        return _pysimple.f90wrap_params__get__trace_time()
    
    @trace_time.setter
    def trace_time(self, trace_time):
        _pysimple.f90wrap_params__set__trace_time(trace_time)
    
    @property
    def ntimstep(self):
        """
        Element ntimstep ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 64
        
        """
        return _pysimple.f90wrap_params__get__ntimstep()
    
    @ntimstep.setter
    def ntimstep(self, ntimstep):
        _pysimple.f90wrap_params__set__ntimstep(ntimstep)
    
    @property
    def npoiper(self):
        """
        Element npoiper ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 64
        
        """
        return _pysimple.f90wrap_params__get__npoiper()
    
    @npoiper.setter
    def npoiper(self, npoiper):
        _pysimple.f90wrap_params__set__npoiper(npoiper)
    
    @property
    def npoiper2(self):
        """
        Element npoiper2 ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 64
        
        """
        return _pysimple.f90wrap_params__get__npoiper2()
    
    @npoiper2.setter
    def npoiper2(self, npoiper2):
        _pysimple.f90wrap_params__set__npoiper2(npoiper2)
    
    @property
    def n_e(self):
        """
        Element n_e ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 64
        
        """
        return _pysimple.f90wrap_params__get__n_e()
    
    @n_e.setter
    def n_e(self, n_e):
        _pysimple.f90wrap_params__set__n_e(n_e)
    
    @property
    def n_d(self):
        """
        Element n_d ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 64
        
        """
        return _pysimple.f90wrap_params__get__n_d()
    
    @n_d.setter
    def n_d(self, n_d):
        _pysimple.f90wrap_params__set__n_d(n_d)
    
    @property
    def v0(self):
        """
        Element v0 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 65
        
        """
        return _pysimple.f90wrap_params__get__v0()
    
    @v0.setter
    def v0(self, v0):
        _pysimple.f90wrap_params__set__v0(v0)
    
    @property
    def nfirstpart(self):
        """
        Element nfirstpart ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 66
        
        """
        return _pysimple.f90wrap_params__get__nfirstpart()
    
    @nfirstpart.setter
    def nfirstpart(self, nfirstpart):
        _pysimple.f90wrap_params__set__nfirstpart(nfirstpart)
    
    @property
    def nlastpart(self):
        """
        Element nlastpart ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 66
        
        """
        return _pysimple.f90wrap_params__get__nlastpart()
    
    @nlastpart.setter
    def nlastpart(self, nlastpart):
        _pysimple.f90wrap_params__set__nlastpart(nlastpart)
    
    @property
    def debug(self):
        """
        Element debug ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 67
        
        """
        return _pysimple.f90wrap_params__get__debug()
    
    @debug.setter
    def debug(self, debug):
        _pysimple.f90wrap_params__set__debug(debug)
    
    @property
    def ierr(self):
        """
        Element ierr ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/params.f90.i \
            line 68
        
        """
        return _pysimple.f90wrap_params__get__ierr()
    
    @ierr.setter
    def ierr(self, ierr):
        _pysimple.f90wrap_params__set__ierr(ierr)
    
    def __str__(self):
        ret = ['<params>{\n']
        ret.append('    npoi : ')
        ret.append(repr(self.npoi))
        ret.append(',\n    l1i : ')
        ret.append(repr(self.l1i))
        ret.append(',\n    nper : ')
        ret.append(repr(self.nper))
        ret.append(',\n    i : ')
        ret.append(repr(self.i))
        ret.append(',\n    ntestpart : ')
        ret.append(repr(self.ntestpart))
        ret.append(',\n    loopskip : ')
        ret.append(repr(self.loopskip))
        ret.append(',\n    iskip : ')
        ret.append(repr(self.iskip))
        ret.append(',\n    dphi : ')
        ret.append(repr(self.dphi))
        ret.append(',\n    phibeg : ')
        ret.append(repr(self.phibeg))
        ret.append(',\n    bmod00 : ')
        ret.append(repr(self.bmod00))
        ret.append(',\n    rlarm : ')
        ret.append(repr(self.rlarm))
        ret.append(',\n    bmax : ')
        ret.append(repr(self.bmax))
        ret.append(',\n    bmin : ')
        ret.append(repr(self.bmin))
        ret.append(',\n    tau : ')
        ret.append(repr(self.tau))
        ret.append(',\n    dtau : ')
        ret.append(repr(self.dtau))
        ret.append(',\n    dtaumin : ')
        ret.append(repr(self.dtaumin))
        ret.append(',\n    xi : ')
        ret.append(repr(self.xi))
        ret.append(',\n    rt0 : ')
        ret.append(repr(self.rt0))
        ret.append(',\n    r0i : ')
        ret.append(repr(self.r0i))
        ret.append(',\n    cbfi : ')
        ret.append(repr(self.cbfi))
        ret.append(',\n    bz0i : ')
        ret.append(repr(self.bz0i))
        ret.append(',\n    bf0 : ')
        ret.append(repr(self.bf0))
        ret.append(',\n    rbig : ')
        ret.append(repr(self.rbig))
        ret.append(',\n    sbeg : ')
        ret.append(repr(self.sbeg))
        ret.append(',\n    thetabeg : ')
        ret.append(repr(self.thetabeg))
        ret.append(',\n    bstart : ')
        ret.append(repr(self.bstart))
        ret.append(',\n    volstart : ')
        ret.append(repr(self.volstart))
        ret.append(',\n    xstart : ')
        ret.append(repr(self.xstart))
        ret.append(',\n    zstart : ')
        ret.append(repr(self.zstart))
        ret.append(',\n    zend : ')
        ret.append(repr(self.zend))
        ret.append(',\n    confpart_trap : ')
        ret.append(repr(self.confpart_trap))
        ret.append(',\n    confpart_pass : ')
        ret.append(repr(self.confpart_pass))
        ret.append(',\n    times_lost : ')
        ret.append(repr(self.times_lost))
        ret.append(',\n    contr_pp : ')
        ret.append(repr(self.contr_pp))
        ret.append(',\n    ibins : ')
        ret.append(repr(self.ibins))
        ret.append(',\n    startmode : ')
        ret.append(repr(self.startmode))
        ret.append(',\n    ntau : ')
        ret.append(repr(self.ntau))
        ret.append(',\n    integmode : ')
        ret.append(repr(self.integmode))
        ret.append(',\n    kpart : ')
        ret.append(repr(self.kpart))
        ret.append(',\n    relerr : ')
        ret.append(repr(self.relerr))
        ret.append(',\n    trap_par : ')
        ret.append(repr(self.trap_par))
        ret.append(',\n    perp_inv : ')
        ret.append(repr(self.perp_inv))
        ret.append(',\n    iclass : ')
        ret.append(repr(self.iclass))
        ret.append(',\n    n_tip_vars : ')
        ret.append(repr(self.n_tip_vars))
        ret.append(',\n    nplagr : ')
        ret.append(repr(self.nplagr))
        ret.append(',\n    nder : ')
        ret.append(repr(self.nder))
        ret.append(',\n    npl_half : ')
        ret.append(repr(self.npl_half))
        ret.append(',\n    norbper : ')
        ret.append(repr(self.norbper))
        ret.append(',\n    nfp : ')
        ret.append(repr(self.nfp))
        ret.append(',\n    fper : ')
        ret.append(repr(self.fper))
        ret.append(',\n    zerolam : ')
        ret.append(repr(self.zerolam))
        ret.append(',\n    tcut : ')
        ret.append(repr(self.tcut))
        ret.append(',\n    ntcut : ')
        ret.append(repr(self.ntcut))
        ret.append(',\n    class_plot : ')
        ret.append(repr(self.class_plot))
        ret.append(',\n    cut_in_per : ')
        ret.append(repr(self.cut_in_per))
        ret.append(',\n    local : ')
        ret.append(repr(self.local))
        ret.append(',\n    fast_class : ')
        ret.append(repr(self.fast_class))
        ret.append(',\n    swcoll : ')
        ret.append(repr(self.swcoll))
        ret.append(',\n    am1 : ')
        ret.append(repr(self.am1))
        ret.append(',\n    am2 : ')
        ret.append(repr(self.am2))
        ret.append(',\n    z1 : ')
        ret.append(repr(self.z1))
        ret.append(',\n    z2 : ')
        ret.append(repr(self.z2))
        ret.append(',\n    densi1 : ')
        ret.append(repr(self.densi1))
        ret.append(',\n    densi2 : ')
        ret.append(repr(self.densi2))
        ret.append(',\n    tempi1 : ')
        ret.append(repr(self.tempi1))
        ret.append(',\n    tempi2 : ')
        ret.append(repr(self.tempi2))
        ret.append(',\n    tempe : ')
        ret.append(repr(self.tempe))
        ret.append(',\n    dchichi : ')
        ret.append(repr(self.dchichi))
        ret.append(',\n    slowrate : ')
        ret.append(repr(self.slowrate))
        ret.append(',\n    dchichi_norm : ')
        ret.append(repr(self.dchichi_norm))
        ret.append(',\n    slowrate_norm : ')
        ret.append(repr(self.slowrate_norm))
        ret.append(',\n    deterministic : ')
        ret.append(repr(self.deterministic))
        ret.append(',\n    notrace_passing : ')
        ret.append(repr(self.notrace_passing))
        ret.append(',\n    face_al : ')
        ret.append(repr(self.face_al))
        ret.append(',\n    trace_time : ')
        ret.append(repr(self.trace_time))
        ret.append(',\n    ntimstep : ')
        ret.append(repr(self.ntimstep))
        ret.append(',\n    npoiper : ')
        ret.append(repr(self.npoiper))
        ret.append(',\n    npoiper2 : ')
        ret.append(repr(self.npoiper2))
        ret.append(',\n    n_e : ')
        ret.append(repr(self.n_e))
        ret.append(',\n    n_d : ')
        ret.append(repr(self.n_d))
        ret.append(',\n    v0 : ')
        ret.append(repr(self.v0))
        ret.append(',\n    nfirstpart : ')
        ret.append(repr(self.nfirstpart))
        ret.append(',\n    nlastpart : ')
        ret.append(repr(self.nlastpart))
        ret.append(',\n    debug : ')
        ret.append(repr(self.debug))
        ret.append(',\n    ierr : ')
        ret.append(repr(self.ierr))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

params = Params()

class Chamb_Mod(f90wrap.runtime.FortranModule):
    """
    Module chamb_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 6-8
    
    """
    @property
    def rnegflag(self):
        """
        Element rnegflag ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 7
        
        """
        return _pysimple.f90wrap_chamb_mod__get__rnegflag()
    
    @rnegflag.setter
    def rnegflag(self, rnegflag):
        _pysimple.f90wrap_chamb_mod__set__rnegflag(rnegflag)
    
    def __str__(self):
        ret = ['<chamb_mod>{\n']
        ret.append('    rnegflag : ')
        ret.append(repr(self.rnegflag))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

chamb_mod = Chamb_Mod()

class Parmot_Mod(f90wrap.runtime.FortranModule):
    """
    Module parmot_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 11-12
    
    """
    @property
    def rmu(self):
        """
        Element rmu ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 12
        
        """
        return _pysimple.f90wrap_parmot_mod__get__rmu()
    
    @rmu.setter
    def rmu(self, rmu):
        _pysimple.f90wrap_parmot_mod__set__rmu(rmu)
    
    @property
    def ro0(self):
        """
        Element ro0 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 12
        
        """
        return _pysimple.f90wrap_parmot_mod__get__ro0()
    
    @ro0.setter
    def ro0(self, ro0):
        _pysimple.f90wrap_parmot_mod__set__ro0(ro0)
    
    @property
    def eeff(self):
        """
        Element eeff ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 12
        
        """
        return _pysimple.f90wrap_parmot_mod__get__eeff()
    
    @eeff.setter
    def eeff(self, eeff):
        _pysimple.f90wrap_parmot_mod__set__eeff(eeff)
    
    def __str__(self):
        ret = ['<parmot_mod>{\n']
        ret.append('    rmu : ')
        ret.append(repr(self.rmu))
        ret.append(',\n    ro0 : ')
        ret.append(repr(self.ro0))
        ret.append(',\n    eeff : ')
        ret.append(repr(self.eeff))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

parmot_mod = Parmot_Mod()

class New_Vmec_Stuff_Mod(f90wrap.runtime.FortranModule):
    """
    Module new_vmec_stuff_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 15-31
    
    """
    @property
    def nsurfm(self):
        """
        Element nsurfm ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 17
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__nsurfm()
    
    @nsurfm.setter
    def nsurfm(self, nsurfm):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__nsurfm(nsurfm)
    
    @property
    def nstrm(self):
        """
        Element nstrm ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 17
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__nstrm()
    
    @nstrm.setter
    def nstrm(self, nstrm):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__nstrm(nstrm)
    
    @property
    def nper(self):
        """
        Element nper ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 17
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__nper()
    
    @nper.setter
    def nper(self, nper):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__nper(nper)
    
    @property
    def kpar(self):
        """
        Element kpar ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 17
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__kpar()
    
    @kpar.setter
    def kpar(self, kpar):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__kpar(kpar)
    
    @property
    def multharm(self):
        """
        Element multharm ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 18
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__multharm()
    
    @multharm.setter
    def multharm(self, multharm):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__multharm(multharm)
    
    @property
    def n_theta(self):
        """
        Element n_theta ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 18
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__n_theta()
    
    @n_theta.setter
    def n_theta(self, n_theta):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__n_theta(n_theta)
    
    @property
    def n_phi(self):
        """
        Element n_phi ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 18
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__n_phi()
    
    @n_phi.setter
    def n_phi(self, n_phi):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__n_phi(n_phi)
    
    @property
    def ns_a(self):
        """
        Element ns_a ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 19
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__ns_a()
    
    @ns_a.setter
    def ns_a(self, ns_a):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__ns_a(ns_a)
    
    @property
    def ns_s(self):
        """
        Element ns_s ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 20
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__ns_s()
    
    @ns_s.setter
    def ns_s(self, ns_s):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__ns_s(ns_s)
    
    @property
    def ns_tp(self):
        """
        Element ns_tp ftype=integer           pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 21
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__ns_tp()
    
    @ns_tp.setter
    def ns_tp(self, ns_tp):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__ns_tp(ns_tp)
    
    @property
    def rmajor(self):
        """
        Element rmajor ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 22
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__rmajor()
    
    @rmajor.setter
    def rmajor(self, rmajor):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__rmajor(rmajor)
    
    @property
    def h_theta(self):
        """
        Element h_theta ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 22
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__h_theta()
    
    @h_theta.setter
    def h_theta(self, h_theta):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__h_theta(h_theta)
    
    @property
    def h_phi(self):
        """
        Element h_phi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 22
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__h_phi()
    
    @h_phi.setter
    def h_phi(self, h_phi):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__h_phi(h_phi)
    
    @property
    def vmec_b_scale(self):
        """
        Element vmec_b_scale ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 23
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__vmec_b_scale()
    
    @vmec_b_scale.setter
    def vmec_b_scale(self, vmec_b_scale):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__vmec_b_scale(vmec_b_scale)
    
    @property
    def vmec_rz_scale(self):
        """
        Element vmec_rz_scale ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 23
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__vmec_rz_scale()
    
    @vmec_rz_scale.setter
    def vmec_rz_scale(self, vmec_rz_scale):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__vmec_rz_scale(vmec_rz_scale)
    
    @property
    def axm(self):
        """
        Element axm ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__axm(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            axm = self._arrays[array_handle]
        else:
            axm = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__axm)
            self._arrays[array_handle] = axm
        return axm
    
    @axm.setter
    def axm(self, axm):
        self.axm[...] = axm
    
    @property
    def axn(self):
        """
        Element axn ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__axn(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            axn = self._arrays[array_handle]
        else:
            axn = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__axn)
            self._arrays[array_handle] = axn
        return axn
    
    @axn.setter
    def axn(self, axn):
        self.axn[...] = axn
    
    @property
    def soa(self):
        """
        Element soa ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 24
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__soa(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            soa = self._arrays[array_handle]
        else:
            soa = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__soa)
            self._arrays[array_handle] = soa
        return soa
    
    @soa.setter
    def soa(self, soa):
        self.soa[...] = soa
    
    @property
    def aiota(self):
        """
        Element aiota ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__aiota(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            aiota = self._arrays[array_handle]
        else:
            aiota = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__aiota)
            self._arrays[array_handle] = aiota
        return aiota
    
    @aiota.setter
    def aiota(self, aiota):
        self.aiota[...] = aiota
    
    @property
    def s(self):
        """
        Element s ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__s(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s = self._arrays[array_handle]
        else:
            s = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__s)
            self._arrays[array_handle] = s
        return s
    
    @s.setter
    def s(self, s):
        self.s[...] = s
    
    @property
    def sps(self):
        """
        Element sps ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__sps(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sps = self._arrays[array_handle]
        else:
            sps = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__sps)
            self._arrays[array_handle] = sps
        return sps
    
    @sps.setter
    def sps(self, sps):
        self.sps[...] = sps
    
    @property
    def phi(self):
        """
        Element phi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__phi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            phi = self._arrays[array_handle]
        else:
            phi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__phi)
            self._arrays[array_handle] = phi
        return phi
    
    @phi.setter
    def phi(self, phi):
        self.phi[...] = phi
    
    @property
    def almnc(self):
        """
        Element almnc ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__almnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            almnc = self._arrays[array_handle]
        else:
            almnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__almnc)
            self._arrays[array_handle] = almnc
        return almnc
    
    @almnc.setter
    def almnc(self, almnc):
        self.almnc[...] = almnc
    
    @property
    def rmnc(self):
        """
        Element rmnc ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__rmnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rmnc = self._arrays[array_handle]
        else:
            rmnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__rmnc)
            self._arrays[array_handle] = rmnc
        return rmnc
    
    @rmnc.setter
    def rmnc(self, rmnc):
        self.rmnc[...] = rmnc
    
    @property
    def zmnc(self):
        """
        Element zmnc ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 26
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__zmnc(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zmnc = self._arrays[array_handle]
        else:
            zmnc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__zmnc)
            self._arrays[array_handle] = zmnc
        return zmnc
    
    @zmnc.setter
    def zmnc(self, zmnc):
        self.zmnc[...] = zmnc
    
    @property
    def almns(self):
        """
        Element almns ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__almns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            almns = self._arrays[array_handle]
        else:
            almns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__almns)
            self._arrays[array_handle] = almns
        return almns
    
    @almns.setter
    def almns(self, almns):
        self.almns[...] = almns
    
    @property
    def rmns(self):
        """
        Element rmns ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__rmns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            rmns = self._arrays[array_handle]
        else:
            rmns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__rmns)
            self._arrays[array_handle] = rmns
        return rmns
    
    @rmns.setter
    def rmns(self, rmns):
        self.rmns[...] = rmns
    
    @property
    def zmns(self):
        """
        Element zmns ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 27
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__zmns(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            zmns = self._arrays[array_handle]
        else:
            zmns = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__zmns)
            self._arrays[array_handle] = zmns
        return zmns
    
    @zmns.setter
    def zmns(self, zmns):
        self.zmns[...] = zmns
    
    @property
    def sr(self):
        """
        Element sr ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__sr(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sr = self._arrays[array_handle]
        else:
            sr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__sr)
            self._arrays[array_handle] = sr
        return sr
    
    @sr.setter
    def sr(self, sr):
        self.sr[...] = sr
    
    @property
    def sz(self):
        """
        Element sz ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__sz(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sz = self._arrays[array_handle]
        else:
            sz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__sz)
            self._arrays[array_handle] = sz
        return sz
    
    @sz.setter
    def sz(self, sz):
        self.sz[...] = sz
    
    @property
    def slam(self):
        """
        Element slam ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 29
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_new_vmec_stuff_mod__array__slam(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            slam = self._arrays[array_handle]
        else:
            slam = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_new_vmec_stuff_mod__array__slam)
            self._arrays[array_handle] = slam
        return slam
    
    @slam.setter
    def slam(self, slam):
        self.slam[...] = slam
    
    @property
    def old_axis_healing(self):
        """
        Element old_axis_healing ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 30
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__old_axis_healing()
    
    @old_axis_healing.setter
    def old_axis_healing(self, old_axis_healing):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__old_axis_healing(old_axis_healing)
    
    @property
    def old_axis_healing_boundary(self):
        """
        Element old_axis_healing_boundary ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 31
        
        """
        return _pysimple.f90wrap_new_vmec_stuff_mod__get__old_axis_healing_boundary()
    
    @old_axis_healing_boundary.setter
    def old_axis_healing_boundary(self, old_axis_healing_boundary):
        _pysimple.f90wrap_new_vmec_stuff_mod__set__old_axis_healing_boundary(old_axis_healing_boundary)
    
    def __str__(self):
        ret = ['<new_vmec_stuff_mod>{\n']
        ret.append('    nsurfm : ')
        ret.append(repr(self.nsurfm))
        ret.append(',\n    nstrm : ')
        ret.append(repr(self.nstrm))
        ret.append(',\n    nper : ')
        ret.append(repr(self.nper))
        ret.append(',\n    kpar : ')
        ret.append(repr(self.kpar))
        ret.append(',\n    multharm : ')
        ret.append(repr(self.multharm))
        ret.append(',\n    n_theta : ')
        ret.append(repr(self.n_theta))
        ret.append(',\n    n_phi : ')
        ret.append(repr(self.n_phi))
        ret.append(',\n    ns_a : ')
        ret.append(repr(self.ns_a))
        ret.append(',\n    ns_s : ')
        ret.append(repr(self.ns_s))
        ret.append(',\n    ns_tp : ')
        ret.append(repr(self.ns_tp))
        ret.append(',\n    rmajor : ')
        ret.append(repr(self.rmajor))
        ret.append(',\n    h_theta : ')
        ret.append(repr(self.h_theta))
        ret.append(',\n    h_phi : ')
        ret.append(repr(self.h_phi))
        ret.append(',\n    vmec_b_scale : ')
        ret.append(repr(self.vmec_b_scale))
        ret.append(',\n    vmec_rz_scale : ')
        ret.append(repr(self.vmec_rz_scale))
        ret.append(',\n    axm : ')
        ret.append(repr(self.axm))
        ret.append(',\n    axn : ')
        ret.append(repr(self.axn))
        ret.append(',\n    soa : ')
        ret.append(repr(self.soa))
        ret.append(',\n    aiota : ')
        ret.append(repr(self.aiota))
        ret.append(',\n    s : ')
        ret.append(repr(self.s))
        ret.append(',\n    sps : ')
        ret.append(repr(self.sps))
        ret.append(',\n    phi : ')
        ret.append(repr(self.phi))
        ret.append(',\n    almnc : ')
        ret.append(repr(self.almnc))
        ret.append(',\n    rmnc : ')
        ret.append(repr(self.rmnc))
        ret.append(',\n    zmnc : ')
        ret.append(repr(self.zmnc))
        ret.append(',\n    almns : ')
        ret.append(repr(self.almns))
        ret.append(',\n    rmns : ')
        ret.append(repr(self.rmns))
        ret.append(',\n    zmns : ')
        ret.append(repr(self.zmns))
        ret.append(',\n    sr : ')
        ret.append(repr(self.sr))
        ret.append(',\n    sz : ')
        ret.append(repr(self.sz))
        ret.append(',\n    slam : ')
        ret.append(repr(self.slam))
        ret.append(',\n    old_axis_healing : ')
        ret.append(repr(self.old_axis_healing))
        ret.append(',\n    old_axis_healing_boundary : ')
        ret.append(repr(self.old_axis_healing_boundary))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

new_vmec_stuff_mod = New_Vmec_Stuff_Mod()

class Vector_Potentail_Mod(f90wrap.runtime.FortranModule):
    """
    Module vector_potentail_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 34-37
    
    """
    @property
    def ns(self):
        """
        Element ns ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 35
        
        """
        return _pysimple.f90wrap_vector_potentail_mod__get__ns()
    
    @ns.setter
    def ns(self, ns):
        _pysimple.f90wrap_vector_potentail_mod__set__ns(ns)
    
    @property
    def hs(self):
        """
        Element hs ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 36
        
        """
        return _pysimple.f90wrap_vector_potentail_mod__get__hs()
    
    @hs.setter
    def hs(self, hs):
        _pysimple.f90wrap_vector_potentail_mod__set__hs(hs)
    
    @property
    def torflux(self):
        """
        Element torflux ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 36
        
        """
        return _pysimple.f90wrap_vector_potentail_mod__get__torflux()
    
    @torflux.setter
    def torflux(self, torflux):
        _pysimple.f90wrap_vector_potentail_mod__set__torflux(torflux)
    
    @property
    def sa_phi(self):
        """
        Element sa_phi ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_vector_potentail_mod__array__sa_phi(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sa_phi = self._arrays[array_handle]
        else:
            sa_phi = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_vector_potentail_mod__array__sa_phi)
            self._arrays[array_handle] = sa_phi
        return sa_phi
    
    @sa_phi.setter
    def sa_phi(self, sa_phi):
        self.sa_phi[...] = sa_phi
    
    def __str__(self):
        ret = ['<vector_potentail_mod>{\n']
        ret.append('    ns : ')
        ret.append(repr(self.ns))
        ret.append(',\n    hs : ')
        ret.append(repr(self.hs))
        ret.append(',\n    torflux : ')
        ret.append(repr(self.torflux))
        ret.append(',\n    sa_phi : ')
        ret.append(repr(self.sa_phi))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

vector_potentail_mod = Vector_Potentail_Mod()

class Canonical_Coordinates_Mod(f90wrap.runtime.FortranModule):
    """
    Module canonical_coordinates_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 40-53
    
    """
    @property
    def ns_max(self):
        """
        Element ns_max ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 41
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__ns_max()
    
    @property
    def n_qua(self):
        """
        Element n_qua ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 41
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__n_qua()
    
    @property
    def ns_s_c(self):
        """
        Element ns_s_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 43
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__ns_s_c()
    
    @ns_s_c.setter
    def ns_s_c(self, ns_s_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__ns_s_c(ns_s_c)
    
    @property
    def ns_tp_c(self):
        """
        Element ns_tp_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 43
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__ns_tp_c()
    
    @ns_tp_c.setter
    def ns_tp_c(self, ns_tp_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__ns_tp_c(ns_tp_c)
    
    @property
    def ns_c(self):
        """
        Element ns_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 44
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__ns_c()
    
    @ns_c.setter
    def ns_c(self, ns_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__ns_c(ns_c)
    
    @property
    def n_theta_c(self):
        """
        Element n_theta_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 44
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__n_theta_c()
    
    @n_theta_c.setter
    def n_theta_c(self, n_theta_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__n_theta_c(n_theta_c)
    
    @property
    def n_phi_c(self):
        """
        Element n_phi_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 44
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__n_phi_c()
    
    @n_phi_c.setter
    def n_phi_c(self, n_phi_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__n_phi_c(n_phi_c)
    
    @property
    def nh_stencil(self):
        """
        Element nh_stencil ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 45
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__nh_stencil()
    
    @property
    def hs_c(self):
        """
        Element hs_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 46
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__hs_c()
    
    @hs_c.setter
    def hs_c(self, hs_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__hs_c(hs_c)
    
    @property
    def h_theta_c(self):
        """
        Element h_theta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 46
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__h_theta_c()
    
    @h_theta_c.setter
    def h_theta_c(self, h_theta_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__h_theta_c(h_theta_c)
    
    @property
    def h_phi_c(self):
        """
        Element h_phi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 46
        
        """
        return _pysimple.f90wrap_canonical_coordinates_mod__get__h_phi_c()
    
    @h_phi_c.setter
    def h_phi_c(self, h_phi_c):
        _pysimple.f90wrap_canonical_coordinates_mod__set__h_phi_c(h_phi_c)
    
    @property
    def g_c(self):
        """
        Element g_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__g_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            g_c = self._arrays[array_handle]
        else:
            g_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__g_c)
            self._arrays[array_handle] = g_c
        return g_c
    
    @g_c.setter
    def g_c(self, g_c):
        self.g_c[...] = g_c
    
    @property
    def sqg_c(self):
        """
        Element sqg_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__sqg_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            sqg_c = self._arrays[array_handle]
        else:
            sqg_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__sqg_c)
            self._arrays[array_handle] = sqg_c
        return sqg_c
    
    @sqg_c.setter
    def sqg_c(self, sqg_c):
        self.sqg_c[...] = sqg_c
    
    @property
    def b_vartheta_c(self):
        """
        Element b_vartheta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__b_vartheta_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            b_vartheta_c = self._arrays[array_handle]
        else:
            b_vartheta_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__b_vartheta_c)
            self._arrays[array_handle] = b_vartheta_c
        return b_vartheta_c
    
    @b_vartheta_c.setter
    def b_vartheta_c(self, b_vartheta_c):
        self.b_vartheta_c[...] = b_vartheta_c
    
    @property
    def b_varphi_c(self):
        """
        Element b_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__b_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            b_varphi_c = self._arrays[array_handle]
        else:
            b_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__b_varphi_c)
            self._arrays[array_handle] = b_varphi_c
        return b_varphi_c
    
    @b_varphi_c.setter
    def b_varphi_c(self, b_varphi_c):
        self.b_varphi_c[...] = b_varphi_c
    
    @property
    def a_vartheta_c(self):
        """
        Element a_vartheta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__a_vartheta_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            a_vartheta_c = self._arrays[array_handle]
        else:
            a_vartheta_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__a_vartheta_c)
            self._arrays[array_handle] = a_vartheta_c
        return a_vartheta_c
    
    @a_vartheta_c.setter
    def a_vartheta_c(self, a_vartheta_c):
        self.a_vartheta_c[...] = a_vartheta_c
    
    @property
    def a_varphi_c(self):
        """
        Element a_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__a_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            a_varphi_c = self._arrays[array_handle]
        else:
            a_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__a_varphi_c)
            self._arrays[array_handle] = a_varphi_c
        return a_varphi_c
    
    @a_varphi_c.setter
    def a_varphi_c(self, a_varphi_c):
        self.a_varphi_c[...] = a_varphi_c
    
    @property
    def delta_varphi_c(self):
        """
        Element delta_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__delta_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            delta_varphi_c = self._arrays[array_handle]
        else:
            delta_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__delta_varphi_c)
            self._arrays[array_handle] = delta_varphi_c
        return delta_varphi_c
    
    @delta_varphi_c.setter
    def delta_varphi_c(self, delta_varphi_c):
        self.delta_varphi_c[...] = delta_varphi_c
    
    @property
    def chi_gauge(self):
        """
        Element chi_gauge ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__chi_gauge(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            chi_gauge = self._arrays[array_handle]
        else:
            chi_gauge = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__chi_gauge)
            self._arrays[array_handle] = chi_gauge
        return chi_gauge
    
    @chi_gauge.setter
    def chi_gauge(self, chi_gauge):
        self.chi_gauge[...] = chi_gauge
    
    @property
    def bmod_c(self):
        """
        Element bmod_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 49
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__bmod_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bmod_c = self._arrays[array_handle]
        else:
            bmod_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__bmod_c)
            self._arrays[array_handle] = bmod_c
        return bmod_c
    
    @bmod_c.setter
    def bmod_c(self, bmod_c):
        self.bmod_c[...] = bmod_c
    
    @property
    def derf1(self):
        """
        Element derf1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__derf1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf1 = self._arrays[array_handle]
        else:
            derf1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__derf1)
            self._arrays[array_handle] = derf1
        return derf1
    
    @derf1.setter
    def derf1(self, derf1):
        self.derf1[...] = derf1
    
    @property
    def derf2(self):
        """
        Element derf2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__derf2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf2 = self._arrays[array_handle]
        else:
            derf2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__derf2)
            self._arrays[array_handle] = derf2
        return derf2
    
    @derf2.setter
    def derf2(self, derf2):
        self.derf2[...] = derf2
    
    @property
    def derf3(self):
        """
        Element derf3 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 51
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__derf3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf3 = self._arrays[array_handle]
        else:
            derf3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__derf3)
            self._arrays[array_handle] = derf3
        return derf3
    
    @derf3.setter
    def derf3(self, derf3):
        self.derf3[...] = derf3
    
    @property
    def s_g_c(self):
        """
        Element s_g_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 52
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__s_g_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_g_c = self._arrays[array_handle]
        else:
            s_g_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__s_g_c)
            self._arrays[array_handle] = s_g_c
        return s_g_c
    
    @s_g_c.setter
    def s_g_c(self, s_g_c):
        self.s_g_c[...] = s_g_c
    
    @property
    def s_sqg_bt_bp(self):
        """
        Element s_sqg_bt_bp ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 53
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_mod__array__s_sqg_bt_bp(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_sqg_bt_bp = self._arrays[array_handle]
        else:
            s_sqg_bt_bp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_mod__array__s_sqg_bt_bp)
            self._arrays[array_handle] = s_sqg_bt_bp
        return s_sqg_bt_bp
    
    @s_sqg_bt_bp.setter
    def s_sqg_bt_bp(self, s_sqg_bt_bp):
        self.s_sqg_bt_bp[...] = s_sqg_bt_bp
    
    def __str__(self):
        ret = ['<canonical_coordinates_mod>{\n']
        ret.append('    ns_max : ')
        ret.append(repr(self.ns_max))
        ret.append(',\n    n_qua : ')
        ret.append(repr(self.n_qua))
        ret.append(',\n    ns_s_c : ')
        ret.append(repr(self.ns_s_c))
        ret.append(',\n    ns_tp_c : ')
        ret.append(repr(self.ns_tp_c))
        ret.append(',\n    ns_c : ')
        ret.append(repr(self.ns_c))
        ret.append(',\n    n_theta_c : ')
        ret.append(repr(self.n_theta_c))
        ret.append(',\n    n_phi_c : ')
        ret.append(repr(self.n_phi_c))
        ret.append(',\n    nh_stencil : ')
        ret.append(repr(self.nh_stencil))
        ret.append(',\n    hs_c : ')
        ret.append(repr(self.hs_c))
        ret.append(',\n    h_theta_c : ')
        ret.append(repr(self.h_theta_c))
        ret.append(',\n    h_phi_c : ')
        ret.append(repr(self.h_phi_c))
        ret.append(',\n    g_c : ')
        ret.append(repr(self.g_c))
        ret.append(',\n    sqg_c : ')
        ret.append(repr(self.sqg_c))
        ret.append(',\n    b_vartheta_c : ')
        ret.append(repr(self.b_vartheta_c))
        ret.append(',\n    b_varphi_c : ')
        ret.append(repr(self.b_varphi_c))
        ret.append(',\n    a_vartheta_c : ')
        ret.append(repr(self.a_vartheta_c))
        ret.append(',\n    a_varphi_c : ')
        ret.append(repr(self.a_varphi_c))
        ret.append(',\n    delta_varphi_c : ')
        ret.append(repr(self.delta_varphi_c))
        ret.append(',\n    chi_gauge : ')
        ret.append(repr(self.chi_gauge))
        ret.append(',\n    bmod_c : ')
        ret.append(repr(self.bmod_c))
        ret.append(',\n    derf1 : ')
        ret.append(repr(self.derf1))
        ret.append(',\n    derf2 : ')
        ret.append(repr(self.derf2))
        ret.append(',\n    derf3 : ')
        ret.append(repr(self.derf3))
        ret.append(',\n    s_g_c : ')
        ret.append(repr(self.s_g_c))
        ret.append(',\n    s_sqg_bt_bp : ')
        ret.append(repr(self.s_sqg_bt_bp))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

canonical_coordinates_mod = Canonical_Coordinates_Mod()

class Canonical_Coordinates_New_Mod(f90wrap.runtime.FortranModule):
    """
    Module canonical_coordinates_new_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 55-66
    
    """
    @property
    def ns_max(self):
        """
        Element ns_max ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 57
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__ns_max()
    
    @property
    def n_qua(self):
        """
        Element n_qua ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 57
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__n_qua()
    
    @property
    def ns_s_c(self):
        """
        Element ns_s_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__ns_s_c()
    
    @ns_s_c.setter
    def ns_s_c(self, ns_s_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__ns_s_c(ns_s_c)
    
    @property
    def ns_tp_c(self):
        """
        Element ns_tp_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 58
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__ns_tp_c()
    
    @ns_tp_c.setter
    def ns_tp_c(self, ns_tp_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__ns_tp_c(ns_tp_c)
    
    @property
    def ns_c(self):
        """
        Element ns_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__ns_c()
    
    @ns_c.setter
    def ns_c(self, ns_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__ns_c(ns_c)
    
    @property
    def n_theta_c(self):
        """
        Element n_theta_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__n_theta_c()
    
    @n_theta_c.setter
    def n_theta_c(self, n_theta_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__n_theta_c(n_theta_c)
    
    @property
    def n_phi_c(self):
        """
        Element n_phi_c ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 59
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__n_phi_c()
    
    @n_phi_c.setter
    def n_phi_c(self, n_phi_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__n_phi_c(n_phi_c)
    
    @property
    def nh_stencil(self):
        """
        Element nh_stencil ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 60
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__nh_stencil()
    
    @property
    def hs_c(self):
        """
        Element hs_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 61
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__hs_c()
    
    @hs_c.setter
    def hs_c(self, hs_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__hs_c(hs_c)
    
    @property
    def h_theta_c(self):
        """
        Element h_theta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 61
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__h_theta_c()
    
    @h_theta_c.setter
    def h_theta_c(self, h_theta_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__h_theta_c(h_theta_c)
    
    @property
    def h_phi_c(self):
        """
        Element h_phi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 61
        
        """
        return _pysimple.f90wrap_canonical_coordinates_new_mod__get__h_phi_c()
    
    @h_phi_c.setter
    def h_phi_c(self, h_phi_c):
        _pysimple.f90wrap_canonical_coordinates_new_mod__set__h_phi_c(h_phi_c)
    
    @property
    def delta_varphi_c(self):
        """
        Element delta_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__delta_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            delta_varphi_c = self._arrays[array_handle]
        else:
            delta_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__delta_varphi_c)
            self._arrays[array_handle] = delta_varphi_c
        return delta_varphi_c
    
    @delta_varphi_c.setter
    def delta_varphi_c(self, delta_varphi_c):
        self.delta_varphi_c[...] = delta_varphi_c
    
    @property
    def h_vartheta_c(self):
        """
        Element h_vartheta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__h_vartheta_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            h_vartheta_c = self._arrays[array_handle]
        else:
            h_vartheta_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__h_vartheta_c)
            self._arrays[array_handle] = h_vartheta_c
        return h_vartheta_c
    
    @h_vartheta_c.setter
    def h_vartheta_c(self, h_vartheta_c):
        self.h_vartheta_c[...] = h_vartheta_c
    
    @property
    def h_varphi_c(self):
        """
        Element h_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__h_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            h_varphi_c = self._arrays[array_handle]
        else:
            h_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__h_varphi_c)
            self._arrays[array_handle] = h_varphi_c
        return h_varphi_c
    
    @h_varphi_c.setter
    def h_varphi_c(self, h_varphi_c):
        self.h_varphi_c[...] = h_varphi_c
    
    @property
    def a_vartheta_c(self):
        """
        Element a_vartheta_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__a_vartheta_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            a_vartheta_c = self._arrays[array_handle]
        else:
            a_vartheta_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__a_vartheta_c)
            self._arrays[array_handle] = a_vartheta_c
        return a_vartheta_c
    
    @a_vartheta_c.setter
    def a_vartheta_c(self, a_vartheta_c):
        self.a_vartheta_c[...] = a_vartheta_c
    
    @property
    def a_varphi_c(self):
        """
        Element a_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__a_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            a_varphi_c = self._arrays[array_handle]
        else:
            a_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__a_varphi_c)
            self._arrays[array_handle] = a_varphi_c
        return a_varphi_c
    
    @a_varphi_c.setter
    def a_varphi_c(self, a_varphi_c):
        self.a_varphi_c[...] = a_varphi_c
    
    @property
    def chi_gauge(self):
        """
        Element chi_gauge ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__chi_gauge(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            chi_gauge = self._arrays[array_handle]
        else:
            chi_gauge = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__chi_gauge)
            self._arrays[array_handle] = chi_gauge
        return chi_gauge
    
    @chi_gauge.setter
    def chi_gauge(self, chi_gauge):
        self.chi_gauge[...] = chi_gauge
    
    @property
    def bmod_c(self):
        """
        Element bmod_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 63
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__bmod_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            bmod_c = self._arrays[array_handle]
        else:
            bmod_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__bmod_c)
            self._arrays[array_handle] = bmod_c
        return bmod_c
    
    @bmod_c.setter
    def bmod_c(self, bmod_c):
        self.bmod_c[...] = bmod_c
    
    @property
    def derf1(self):
        """
        Element derf1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__derf1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf1 = self._arrays[array_handle]
        else:
            derf1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__derf1)
            self._arrays[array_handle] = derf1
        return derf1
    
    @derf1.setter
    def derf1(self, derf1):
        self.derf1[...] = derf1
    
    @property
    def derf2(self):
        """
        Element derf2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__derf2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf2 = self._arrays[array_handle]
        else:
            derf2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__derf2)
            self._arrays[array_handle] = derf2
        return derf2
    
    @derf2.setter
    def derf2(self, derf2):
        self.derf2[...] = derf2
    
    @property
    def derf3(self):
        """
        Element derf3 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 64
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__derf3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf3 = self._arrays[array_handle]
        else:
            derf3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__derf3)
            self._arrays[array_handle] = derf3
        return derf3
    
    @derf3.setter
    def derf3(self, derf3):
        self.derf3[...] = derf3
    
    @property
    def s_delta_varphi_c(self):
        """
        Element s_delta_varphi_c ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 65
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__s_delta_varphi_c(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_delta_varphi_c = self._arrays[array_handle]
        else:
            s_delta_varphi_c = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__s_delta_varphi_c)
            self._arrays[array_handle] = s_delta_varphi_c
        return s_delta_varphi_c
    
    @s_delta_varphi_c.setter
    def s_delta_varphi_c(self, s_delta_varphi_c):
        self.s_delta_varphi_c[...] = s_delta_varphi_c
    
    @property
    def s_bmod_b_a(self):
        """
        Element s_bmod_b_a ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 66
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_canonical_coordinates_new_mod__array__s_bmod_b_a(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_bmod_b_a = self._arrays[array_handle]
        else:
            s_bmod_b_a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_canonical_coordinates_new_mod__array__s_bmod_b_a)
            self._arrays[array_handle] = s_bmod_b_a
        return s_bmod_b_a
    
    @s_bmod_b_a.setter
    def s_bmod_b_a(self, s_bmod_b_a):
        self.s_bmod_b_a[...] = s_bmod_b_a
    
    def __str__(self):
        ret = ['<canonical_coordinates_new_mod>{\n']
        ret.append('    ns_max : ')
        ret.append(repr(self.ns_max))
        ret.append(',\n    n_qua : ')
        ret.append(repr(self.n_qua))
        ret.append(',\n    ns_s_c : ')
        ret.append(repr(self.ns_s_c))
        ret.append(',\n    ns_tp_c : ')
        ret.append(repr(self.ns_tp_c))
        ret.append(',\n    ns_c : ')
        ret.append(repr(self.ns_c))
        ret.append(',\n    n_theta_c : ')
        ret.append(repr(self.n_theta_c))
        ret.append(',\n    n_phi_c : ')
        ret.append(repr(self.n_phi_c))
        ret.append(',\n    nh_stencil : ')
        ret.append(repr(self.nh_stencil))
        ret.append(',\n    hs_c : ')
        ret.append(repr(self.hs_c))
        ret.append(',\n    h_theta_c : ')
        ret.append(repr(self.h_theta_c))
        ret.append(',\n    h_phi_c : ')
        ret.append(repr(self.h_phi_c))
        ret.append(',\n    delta_varphi_c : ')
        ret.append(repr(self.delta_varphi_c))
        ret.append(',\n    h_vartheta_c : ')
        ret.append(repr(self.h_vartheta_c))
        ret.append(',\n    h_varphi_c : ')
        ret.append(repr(self.h_varphi_c))
        ret.append(',\n    a_vartheta_c : ')
        ret.append(repr(self.a_vartheta_c))
        ret.append(',\n    a_varphi_c : ')
        ret.append(repr(self.a_varphi_c))
        ret.append(',\n    chi_gauge : ')
        ret.append(repr(self.chi_gauge))
        ret.append(',\n    bmod_c : ')
        ret.append(repr(self.bmod_c))
        ret.append(',\n    derf1 : ')
        ret.append(repr(self.derf1))
        ret.append(',\n    derf2 : ')
        ret.append(repr(self.derf2))
        ret.append(',\n    derf3 : ')
        ret.append(repr(self.derf3))
        ret.append(',\n    s_delta_varphi_c : ')
        ret.append(repr(self.s_delta_varphi_c))
        ret.append(',\n    s_bmod_b_a : ')
        ret.append(repr(self.s_bmod_b_a))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

canonical_coordinates_new_mod = Canonical_Coordinates_New_Mod()

class Boozer_Coordinates_Mod(f90wrap.runtime.FortranModule):
    """
    Module boozer_coordinates_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 69-89
    
    """
    @property
    def ns_max(self):
        """
        Element ns_max ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 70
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__ns_max()
    
    @property
    def n_qua(self):
        """
        Element n_qua ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 70
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__n_qua()
    
    @property
    def use_b_r(self):
        """
        Element use_b_r ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 72
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__use_b_r()
    
    @use_b_r.setter
    def use_b_r(self, use_b_r):
        _pysimple.f90wrap_boozer_coordinates_mod__set__use_b_r(use_b_r)
    
    @property
    def use_del_tp_b(self):
        """
        Element use_del_tp_b ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 72
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__use_del_tp_b()
    
    @use_del_tp_b.setter
    def use_del_tp_b(self, use_del_tp_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__use_del_tp_b(use_del_tp_b)
    
    @property
    def ns_s_b(self):
        """
        Element ns_s_b ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 74
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__ns_s_b()
    
    @ns_s_b.setter
    def ns_s_b(self, ns_s_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__ns_s_b(ns_s_b)
    
    @property
    def ns_tp_b(self):
        """
        Element ns_tp_b ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 74
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__ns_tp_b()
    
    @ns_tp_b.setter
    def ns_tp_b(self, ns_tp_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__ns_tp_b(ns_tp_b)
    
    @property
    def ns_b(self):
        """
        Element ns_b ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 75
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__ns_b()
    
    @ns_b.setter
    def ns_b(self, ns_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__ns_b(ns_b)
    
    @property
    def n_theta_b(self):
        """
        Element n_theta_b ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 75
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__n_theta_b()
    
    @n_theta_b.setter
    def n_theta_b(self, n_theta_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__n_theta_b(n_theta_b)
    
    @property
    def n_phi_b(self):
        """
        Element n_phi_b ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 75
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__n_phi_b()
    
    @n_phi_b.setter
    def n_phi_b(self, n_phi_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__n_phi_b(n_phi_b)
    
    @property
    def hs_b(self):
        """
        Element hs_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 76
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__hs_b()
    
    @hs_b.setter
    def hs_b(self, hs_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__hs_b(hs_b)
    
    @property
    def h_theta_b(self):
        """
        Element h_theta_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 76
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__h_theta_b()
    
    @h_theta_b.setter
    def h_theta_b(self, h_theta_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__h_theta_b(h_theta_b)
    
    @property
    def h_phi_b(self):
        """
        Element h_phi_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 76
        
        """
        return _pysimple.f90wrap_boozer_coordinates_mod__get__h_phi_b()
    
    @h_phi_b.setter
    def h_phi_b(self, h_phi_b):
        _pysimple.f90wrap_boozer_coordinates_mod__set__h_phi_b(h_phi_b)
    
    @property
    def derf1(self):
        """
        Element derf1 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__derf1(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf1 = self._arrays[array_handle]
        else:
            derf1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__derf1)
            self._arrays[array_handle] = derf1
        return derf1
    
    @derf1.setter
    def derf1(self, derf1):
        self.derf1[...] = derf1
    
    @property
    def derf2(self):
        """
        Element derf2 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__derf2(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf2 = self._arrays[array_handle]
        else:
            derf2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__derf2)
            self._arrays[array_handle] = derf2
        return derf2
    
    @derf2.setter
    def derf2(self, derf2):
        self.derf2[...] = derf2
    
    @property
    def derf3(self):
        """
        Element derf3 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__derf3(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            derf3 = self._arrays[array_handle]
        else:
            derf3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__derf3)
            self._arrays[array_handle] = derf3
        return derf3
    
    @derf3.setter
    def derf3(self, derf3):
        self.derf3[...] = derf3
    
    @property
    def s_bcovar_tp_b(self):
        """
        Element s_bcovar_tp_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 81
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__s_bcovar_tp_b(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_bcovar_tp_b = self._arrays[array_handle]
        else:
            s_bcovar_tp_b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__s_bcovar_tp_b)
            self._arrays[array_handle] = s_bcovar_tp_b
        return s_bcovar_tp_b
    
    @s_bcovar_tp_b.setter
    def s_bcovar_tp_b(self, s_bcovar_tp_b):
        self.s_bcovar_tp_b[...] = s_bcovar_tp_b
    
    @property
    def s_bmod_b(self):
        """
        Element s_bmod_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 84
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__s_bmod_b(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_bmod_b = self._arrays[array_handle]
        else:
            s_bmod_b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__s_bmod_b)
            self._arrays[array_handle] = s_bmod_b
        return s_bmod_b
    
    @s_bmod_b.setter
    def s_bmod_b(self, s_bmod_b):
        self.s_bmod_b[...] = s_bmod_b
    
    @property
    def s_bcovar_r_b(self):
        """
        Element s_bcovar_r_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 84
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__s_bcovar_r_b(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_bcovar_r_b = self._arrays[array_handle]
        else:
            s_bcovar_r_b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__s_bcovar_r_b)
            self._arrays[array_handle] = s_bcovar_r_b
        return s_bcovar_r_b
    
    @s_bcovar_r_b.setter
    def s_bcovar_r_b(self, s_bcovar_r_b):
        self.s_bcovar_r_b[...] = s_bcovar_r_b
    
    @property
    def s_delt_delp_v(self):
        """
        Element s_delt_delp_v ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__s_delt_delp_v(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_delt_delp_v = self._arrays[array_handle]
        else:
            s_delt_delp_v = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__s_delt_delp_v)
            self._arrays[array_handle] = s_delt_delp_v
        return s_delt_delp_v
    
    @s_delt_delp_v.setter
    def s_delt_delp_v(self, s_delt_delp_v):
        self.s_delt_delp_v[...] = s_delt_delp_v
    
    @property
    def s_delt_delp_b(self):
        """
        Element s_delt_delp_b ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_boozer_coordinates_mod__array__s_delt_delp_b(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            s_delt_delp_b = self._arrays[array_handle]
        else:
            s_delt_delp_b = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_boozer_coordinates_mod__array__s_delt_delp_b)
            self._arrays[array_handle] = s_delt_delp_b
        return s_delt_delp_b
    
    @s_delt_delp_b.setter
    def s_delt_delp_b(self, s_delt_delp_b):
        self.s_delt_delp_b[...] = s_delt_delp_b
    
    def __str__(self):
        ret = ['<boozer_coordinates_mod>{\n']
        ret.append('    ns_max : ')
        ret.append(repr(self.ns_max))
        ret.append(',\n    n_qua : ')
        ret.append(repr(self.n_qua))
        ret.append(',\n    use_b_r : ')
        ret.append(repr(self.use_b_r))
        ret.append(',\n    use_del_tp_b : ')
        ret.append(repr(self.use_del_tp_b))
        ret.append(',\n    ns_s_b : ')
        ret.append(repr(self.ns_s_b))
        ret.append(',\n    ns_tp_b : ')
        ret.append(repr(self.ns_tp_b))
        ret.append(',\n    ns_b : ')
        ret.append(repr(self.ns_b))
        ret.append(',\n    n_theta_b : ')
        ret.append(repr(self.n_theta_b))
        ret.append(',\n    n_phi_b : ')
        ret.append(repr(self.n_phi_b))
        ret.append(',\n    hs_b : ')
        ret.append(repr(self.hs_b))
        ret.append(',\n    h_theta_b : ')
        ret.append(repr(self.h_theta_b))
        ret.append(',\n    h_phi_b : ')
        ret.append(repr(self.h_phi_b))
        ret.append(',\n    derf1 : ')
        ret.append(repr(self.derf1))
        ret.append(',\n    derf2 : ')
        ret.append(repr(self.derf2))
        ret.append(',\n    derf3 : ')
        ret.append(repr(self.derf3))
        ret.append(',\n    s_bcovar_tp_b : ')
        ret.append(repr(self.s_bcovar_tp_b))
        ret.append(',\n    s_bmod_b : ')
        ret.append(repr(self.s_bmod_b))
        ret.append(',\n    s_bcovar_r_b : ')
        ret.append(repr(self.s_bcovar_r_b))
        ret.append(',\n    s_delt_delp_v : ')
        ret.append(repr(self.s_delt_delp_v))
        ret.append(',\n    s_delt_delp_b : ')
        ret.append(repr(self.s_delt_delp_b))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

boozer_coordinates_mod = Boozer_Coordinates_Mod()

class Velo_Mod(f90wrap.runtime.FortranModule):
    """
    Module velo_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 92-93
    
    """
    @property
    def isw_field_type(self):
        """
        Element isw_field_type ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 93
        
        """
        return _pysimple.f90wrap_velo_mod__get__isw_field_type()
    
    @isw_field_type.setter
    def isw_field_type(self, isw_field_type):
        _pysimple.f90wrap_velo_mod__set__isw_field_type(isw_field_type)
    
    def __str__(self):
        ret = ['<velo_mod>{\n']
        ret.append('    isw_field_type : ')
        ret.append(repr(self.isw_field_type))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

velo_mod = Velo_Mod()

class Gbpi_Mod(f90wrap.runtime.FortranModule):
    """
    Module gbpi_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 96-98
    
    """
    @property
    def ierrfield(self):
        """
        Element ierrfield ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 98
        
        """
        return _pysimple.f90wrap_gbpi_mod__get__ierrfield()
    
    @ierrfield.setter
    def ierrfield(self, ierrfield):
        _pysimple.f90wrap_gbpi_mod__set__ierrfield(ierrfield)
    
    def __str__(self):
        ret = ['<gbpi_mod>{\n']
        ret.append('    ierrfield : ')
        ret.append(repr(self.ierrfield))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

gbpi_mod = Gbpi_Mod()

class Diag_Mod(f90wrap.runtime.FortranModule):
    """
    Module diag_mod
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
        lines 101-102
    
    """
    @property
    def dodiag(self):
        """
        Element dodiag ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 102
        
        """
        return _pysimple.f90wrap_diag_mod__get__dodiag()
    
    @dodiag.setter
    def dodiag(self, dodiag):
        _pysimple.f90wrap_diag_mod__set__dodiag(dodiag)
    
    @property
    def icounter(self):
        """
        Element icounter ftype=integer(8) pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/canonical_coordinates_mod.f90.i \
            line 103
        
        """
        return _pysimple.f90wrap_diag_mod__get__icounter()
    
    @icounter.setter
    def icounter(self, icounter):
        _pysimple.f90wrap_diag_mod__set__icounter(icounter)
    
    def __str__(self):
        ret = ['<diag_mod>{\n']
        ret.append('    dodiag : ')
        ret.append(repr(self.dodiag))
        ret.append(',\n    icounter : ')
        ret.append(repr(self.icounter))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

diag_mod = Diag_Mod()

class Simple_Bench(f90wrap.runtime.FortranModule):
    """
    Module simple_bench
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
        lines 5-256
    
    """
    @staticmethod
    def init_bench():
        """
        init_bench()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 41-75
        
        
        """
        _pysimple.f90wrap_init_bench()
    
    @staticmethod
    def do_bench():
        """
        do_bench()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 77-116
        
        
        """
        _pysimple.f90wrap_do_bench()
    
    @staticmethod
    def cleanup_bench():
        """
        cleanup_bench()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 118-120
        
        
        """
        _pysimple.f90wrap_cleanup_bench()
    
    @staticmethod
    def test_cuts(nplagr):
        """
        test_cuts(nplagr)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 122-189
        
        Parameters
        ----------
        nplagr : int
        
        """
        _pysimple.f90wrap_test_cuts(nplagr=nplagr)
    
    @staticmethod
    def test_orbit():
        """
        test_orbit()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 191-192
        
        
        """
        _pysimple.f90wrap_test_orbit()
    
    @staticmethod
    def test_single():
        """
        test_single()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 194-202
        
        
        """
        _pysimple.f90wrap_test_single()
    
    @staticmethod
    def test_quasi():
        """
        test_quasi()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 204-214
        
        
        """
        _pysimple.f90wrap_test_quasi()
    
    @staticmethod
    def test_multi():
        """
        test_multi()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 216-224
        
        
        """
        _pysimple.f90wrap_test_multi()
    
    @staticmethod
    def test_multi_quasi():
        """
        test_multi_quasi()
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 226-236
        
        
        """
        _pysimple.f90wrap_test_multi_quasi()
    
    @staticmethod
    def minsqdist(za, zref, result):
        """
        minsqdist(za, zref, result)
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            lines 238-256
        
        Parameters
        ----------
        za : unknown array
        zref : unknown array
        result : unknown array
        
        """
        _pysimple.f90wrap_minsqdist(za=za, zref=zref, result=result)
    
    @property
    def ierr(self):
        """
        Element ierr ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 16
        
        """
        return _pysimple.f90wrap_simple_bench__get__ierr()
    
    @ierr.setter
    def ierr(self, ierr):
        _pysimple.f90wrap_simple_bench__set__ierr(ierr)
    
    @property
    def kt(self):
        """
        Element kt ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 16
        
        """
        return _pysimple.f90wrap_simple_bench__get__kt()
    
    @kt.setter
    def kt(self, kt):
        _pysimple.f90wrap_simple_bench__set__kt(kt)
    
    @property
    def z0(self):
        """
        Element z0 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 17
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_simple_bench__array__z0(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            z0 = self._arrays[array_handle]
        else:
            z0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_simple_bench__array__z0)
            self._arrays[array_handle] = z0
        return z0
    
    @z0.setter
    def z0(self, z0):
        self.z0[...] = z0
    
    @property
    def vpar0(self):
        """
        Element vpar0 ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 17
        
        """
        return _pysimple.f90wrap_simple_bench__get__vpar0()
    
    @vpar0.setter
    def vpar0(self, vpar0):
        _pysimple.f90wrap_simple_bench__set__vpar0(vpar0)
    
    @property
    def dt(self):
        """
        Element dt ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 17
        
        """
        return _pysimple.f90wrap_simple_bench__get__dt()
    
    @dt.setter
    def dt(self, dt):
        _pysimple.f90wrap_simple_bench__set__dt(dt)
    
    @property
    def npoiper2(self):
        """
        Element npoiper2 ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 22
        
        """
        return _pysimple.f90wrap_simple_bench__get__npoiper2()
    
    @npoiper2.setter
    def npoiper2(self, npoiper2):
        _pysimple.f90wrap_simple_bench__set__npoiper2(npoiper2)
    
    @property
    def rbig(self):
        """
        Element rbig ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 23
        
        """
        return _pysimple.f90wrap_simple_bench__get__rbig()
    
    @rbig.setter
    def rbig(self, rbig):
        _pysimple.f90wrap_simple_bench__set__rbig(rbig)
    
    @property
    def dtau(self):
        """
        Element dtau ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 23
        
        """
        return _pysimple.f90wrap_simple_bench__get__dtau()
    
    @dtau.setter
    def dtau(self, dtau):
        _pysimple.f90wrap_simple_bench__set__dtau(dtau)
    
    @property
    def dtaumax(self):
        """
        Element dtaumax ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 23
        
        """
        return _pysimple.f90wrap_simple_bench__get__dtaumax()
    
    @dtaumax.setter
    def dtaumax(self, dtaumax):
        _pysimple.f90wrap_simple_bench__set__dtaumax(dtaumax)
    
    @property
    def nt(self):
        """
        Element nt ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 24
        
        """
        return _pysimple.f90wrap_simple_bench__get__nt()
    
    @nt.setter
    def nt(self, nt):
        _pysimple.f90wrap_simple_bench__set__nt(nt)
    
    @property
    def out(self):
        """
        Element out ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 25
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_simple_bench__array__out(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            out = self._arrays[array_handle]
        else:
            out = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_simple_bench__array__out)
            self._arrays[array_handle] = out
        return out
    
    @out.setter
    def out(self, out):
        self.out[...] = out
    
    @property
    def starttime(self):
        """
        Element starttime ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 26
        
        """
        return _pysimple.f90wrap_simple_bench__get__starttime()
    
    @starttime.setter
    def starttime(self, starttime):
        _pysimple.f90wrap_simple_bench__set__starttime(starttime)
    
    @property
    def endtime(self):
        """
        Element endtime ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 26
        
        """
        return _pysimple.f90wrap_simple_bench__get__endtime()
    
    @endtime.setter
    def endtime(self, endtime):
        _pysimple.f90wrap_simple_bench__set__endtime(endtime)
    
    @property
    def multi(self):
        """
        Element multi ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 27
        
        """
        return _pysimple.f90wrap_simple_bench__get__multi()
    
    @multi.setter
    def multi(self, multi):
        _pysimple.f90wrap_simple_bench__set__multi(multi)
    
    @property
    def quasi(self):
        """
        Element quasi ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 28
        
        """
        return _pysimple.f90wrap_simple_bench__get__quasi()
    
    @quasi.setter
    def quasi(self, quasi):
        _pysimple.f90wrap_simple_bench__set__quasi(quasi)
    
    @property
    def tok(self):
        """
        Element tok ftype=logical pytype=bool
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 29
        
        """
        return _pysimple.f90wrap_simple_bench__get__tok()
    
    @tok.setter
    def tok(self, tok):
        _pysimple.f90wrap_simple_bench__set__tok(tok)
    
    @property
    def integ_mode(self):
        """
        Element integ_mode ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 30
        
        """
        return _pysimple.f90wrap_simple_bench__get__integ_mode()
    
    @integ_mode.setter
    def integ_mode(self, integ_mode):
        _pysimple.f90wrap_simple_bench__set__integ_mode(integ_mode)
    
    @property
    def nlag(self):
        """
        Element nlag ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 31
        
        """
        return _pysimple.f90wrap_simple_bench__get__nlag()
    
    @nlag.setter
    def nlag(self, nlag):
        _pysimple.f90wrap_simple_bench__set__nlag(nlag)
    
    @property
    def nplagr_invar(self):
        """
        Element nplagr_invar ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 32
        
        """
        return _pysimple.f90wrap_simple_bench__get__nplagr_invar()
    
    @nplagr_invar.setter
    def nplagr_invar(self, nplagr_invar):
        _pysimple.f90wrap_simple_bench__set__nplagr_invar(nplagr_invar)
    
    @property
    def ncut(self):
        """
        Element ncut ftype=integer  pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 33
        
        """
        return _pysimple.f90wrap_simple_bench__get__ncut()
    
    @ncut.setter
    def ncut(self, ncut):
        _pysimple.f90wrap_simple_bench__set__ncut(ncut)
    
    @property
    def infile(self):
        """
        Element infile ftype=character(len=*) pytype=str
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 34
        
        """
        return _pysimple.f90wrap_simple_bench__get__infile()
    
    @property
    def outfile(self):
        """
        Element outfile ftype=character(len=*) pytype=str
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 35
        
        """
        return _pysimple.f90wrap_simple_bench__get__outfile()
    
    @property
    def n_tip_vars(self):
        """
        Element n_tip_vars ftype=integer pytype=int
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 36
        
        """
        return _pysimple.f90wrap_simple_bench__get__n_tip_vars()
    
    @property
    def var_cut(self):
        """
        Element var_cut ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 37
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _pysimple.f90wrap_simple_bench__array__var_cut(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            var_cut = self._arrays[array_handle]
        else:
            var_cut = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _pysimple.f90wrap_simple_bench__array__var_cut)
            self._arrays[array_handle] = var_cut
        return var_cut
    
    @var_cut.setter
    def var_cut(self, var_cut):
        self.var_cut[...] = var_cut
    
    @property
    def taub(self):
        """
        Element taub ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 38
        
        """
        return _pysimple.f90wrap_simple_bench__get__taub()
    
    @taub.setter
    def taub(self, taub):
        _pysimple.f90wrap_simple_bench__set__taub(taub)
    
    @property
    def rtol(self):
        """
        Element rtol ftype=double precision pytype=unknown
        
        
        Defined at \
            /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/bench.f90.i \
            line 39
        
        """
        return _pysimple.f90wrap_simple_bench__get__rtol()
    
    @rtol.setter
    def rtol(self, rtol):
        _pysimple.f90wrap_simple_bench__set__rtol(rtol)
    
    def __str__(self):
        ret = ['<simple_bench>{\n']
        ret.append('    ierr : ')
        ret.append(repr(self.ierr))
        ret.append(',\n    kt : ')
        ret.append(repr(self.kt))
        ret.append(',\n    z0 : ')
        ret.append(repr(self.z0))
        ret.append(',\n    vpar0 : ')
        ret.append(repr(self.vpar0))
        ret.append(',\n    dt : ')
        ret.append(repr(self.dt))
        ret.append(',\n    npoiper2 : ')
        ret.append(repr(self.npoiper2))
        ret.append(',\n    rbig : ')
        ret.append(repr(self.rbig))
        ret.append(',\n    dtau : ')
        ret.append(repr(self.dtau))
        ret.append(',\n    dtaumax : ')
        ret.append(repr(self.dtaumax))
        ret.append(',\n    nt : ')
        ret.append(repr(self.nt))
        ret.append(',\n    out : ')
        ret.append(repr(self.out))
        ret.append(',\n    starttime : ')
        ret.append(repr(self.starttime))
        ret.append(',\n    endtime : ')
        ret.append(repr(self.endtime))
        ret.append(',\n    multi : ')
        ret.append(repr(self.multi))
        ret.append(',\n    quasi : ')
        ret.append(repr(self.quasi))
        ret.append(',\n    tok : ')
        ret.append(repr(self.tok))
        ret.append(',\n    integ_mode : ')
        ret.append(repr(self.integ_mode))
        ret.append(',\n    nlag : ')
        ret.append(repr(self.nlag))
        ret.append(',\n    nplagr_invar : ')
        ret.append(repr(self.nplagr_invar))
        ret.append(',\n    ncut : ')
        ret.append(repr(self.ncut))
        ret.append(',\n    infile : ')
        ret.append(repr(self.infile))
        ret.append(',\n    outfile : ')
        ret.append(repr(self.outfile))
        ret.append(',\n    n_tip_vars : ')
        ret.append(repr(self.n_tip_vars))
        ret.append(',\n    var_cut : ')
        ret.append(repr(self.var_cut))
        ret.append(',\n    taub : ')
        ret.append(repr(self.taub))
        ret.append(',\n    rtol : ')
        ret.append(repr(self.rtol))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

simple_bench = Simple_Bench()

def vmec_to_cyl(s, theta, varphi):
    """
    rcyl, zcyl = vmec_to_cyl(s, theta, varphi)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/simple.f90.i \
        lines 1003-1011
    
    Parameters
    ----------
    s : unknown
    theta : unknown
    varphi : unknown
    
    Returns
    -------
    rcyl : unknown
    zcyl : unknown
    
    """
    rcyl, zcyl = _pysimple.f90wrap_vmec_to_cyl(s=s, theta=theta, varphi=varphi)
    return rcyl, zcyl

def get_canonical_coordinates():
    """
    get_canonical_coordinates()
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 17-220
    
    
    """
    _pysimple.f90wrap_get_canonical_coordinates()

def rhs_cancoord(r, y, dy):
    """
    rhs_cancoord(r, y, dy)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 225-271
    
    Parameters
    ----------
    r : unknown
    y : unknown array
    dy : unknown array
    
    """
    _pysimple.f90wrap_rhs_cancoord(r=r, y=y, dy=dy)

def spline_can_coord(fullset):
    """
    spline_can_coord(fullset)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 276-414
    
    Parameters
    ----------
    fullset : bool
    
    """
    _pysimple.f90wrap_spline_can_coord(fullset=fullset)

def splint_can_coord(fullset, mode_secders, r, vartheta_c, varphi_c, a_theta, \
    a_phi, da_theta_dr, da_phi_dr, d2a_phi_dr2, d3a_phi_dr3, sqg_c, dsqg_c_dr, \
    dsqg_c_dt, dsqg_c_dp, b_vartheta_c, db_vartheta_c_dr, db_vartheta_c_dt, \
    db_vartheta_c_dp, b_varphi_c, db_varphi_c_dr, db_varphi_c_dt, \
    db_varphi_c_dp, d2sqg_rr, d2sqg_rt, d2sqg_rp, d2sqg_tt, d2sqg_tp, d2sqg_pp, \
    d2bth_rr, d2bth_rt, d2bth_rp, d2bth_tt, d2bth_tp, d2bth_pp, d2bph_rr, \
    d2bph_rt, d2bph_rp, d2bph_tt, d2bph_tp, d2bph_pp, g_c):
    """
    splint_can_coord(fullset, mode_secders, r, vartheta_c, varphi_c, a_theta, a_phi, \
        da_theta_dr, da_phi_dr, d2a_phi_dr2, d3a_phi_dr3, sqg_c, dsqg_c_dr, \
        dsqg_c_dt, dsqg_c_dp, b_vartheta_c, db_vartheta_c_dr, db_vartheta_c_dt, \
        db_vartheta_c_dp, b_varphi_c, db_varphi_c_dr, db_varphi_c_dt, \
        db_varphi_c_dp, d2sqg_rr, d2sqg_rt, d2sqg_rp, d2sqg_tt, d2sqg_tp, d2sqg_pp, \
        d2bth_rr, d2bth_rt, d2bth_rp, d2bth_tt, d2bth_tp, d2bth_pp, d2bph_rr, \
        d2bph_rt, d2bph_rp, d2bph_tt, d2bph_tp, d2bph_pp, g_c)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 426-880
    
    Parameters
    ----------
    fullset : bool
    mode_secders : int
    r : unknown
    vartheta_c : unknown
    varphi_c : unknown
    a_theta : unknown
    a_phi : unknown
    da_theta_dr : unknown
    da_phi_dr : unknown
    d2a_phi_dr2 : unknown
    d3a_phi_dr3 : unknown
    sqg_c : unknown
    dsqg_c_dr : unknown
    dsqg_c_dt : unknown
    dsqg_c_dp : unknown
    b_vartheta_c : unknown
    db_vartheta_c_dr : unknown
    db_vartheta_c_dt : unknown
    db_vartheta_c_dp : unknown
    b_varphi_c : unknown
    db_varphi_c_dr : unknown
    db_varphi_c_dt : unknown
    db_varphi_c_dp : unknown
    d2sqg_rr : unknown
    d2sqg_rt : unknown
    d2sqg_rp : unknown
    d2sqg_tt : unknown
    d2sqg_tp : unknown
    d2sqg_pp : unknown
    d2bth_rr : unknown
    d2bth_rt : unknown
    d2bth_rp : unknown
    d2bth_tt : unknown
    d2bth_tp : unknown
    d2bth_pp : unknown
    d2bph_rr : unknown
    d2bph_rt : unknown
    d2bph_rp : unknown
    d2bph_tt : unknown
    d2bph_tp : unknown
    d2bph_pp : unknown
    g_c : unknown
    
    --------------------------------
    -------------------------------
    """
    _pysimple.f90wrap_splint_can_coord(fullset=fullset, mode_secders=mode_secders, \
        r=r, vartheta_c=vartheta_c, varphi_c=varphi_c, a_theta=a_theta, a_phi=a_phi, \
        da_theta_dr=da_theta_dr, da_phi_dr=da_phi_dr, d2a_phi_dr2=d2a_phi_dr2, \
        d3a_phi_dr3=d3a_phi_dr3, sqg_c=sqg_c, dsqg_c_dr=dsqg_c_dr, \
        dsqg_c_dt=dsqg_c_dt, dsqg_c_dp=dsqg_c_dp, b_vartheta_c=b_vartheta_c, \
        db_vartheta_c_dr=db_vartheta_c_dr, db_vartheta_c_dt=db_vartheta_c_dt, \
        db_vartheta_c_dp=db_vartheta_c_dp, b_varphi_c=b_varphi_c, \
        db_varphi_c_dr=db_varphi_c_dr, db_varphi_c_dt=db_varphi_c_dt, \
        db_varphi_c_dp=db_varphi_c_dp, d2sqg_rr=d2sqg_rr, d2sqg_rt=d2sqg_rt, \
        d2sqg_rp=d2sqg_rp, d2sqg_tt=d2sqg_tt, d2sqg_tp=d2sqg_tp, d2sqg_pp=d2sqg_pp, \
        d2bth_rr=d2bth_rr, d2bth_rt=d2bth_rt, d2bth_rp=d2bth_rp, d2bth_tt=d2bth_tt, \
        d2bth_tp=d2bth_tp, d2bth_pp=d2bth_pp, d2bph_rr=d2bph_rr, d2bph_rt=d2bph_rt, \
        d2bph_rp=d2bph_rp, d2bph_tt=d2bph_tt, d2bph_tp=d2bph_tp, d2bph_pp=d2bph_pp, \
        g_c=g_c)

def can_to_vmec(r, vartheta_c_in, varphi_c_in):
    """
    theta_vmec, varphi_vmec = can_to_vmec(r, vartheta_c_in, varphi_c_in)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 885-925
    
    Parameters
    ----------
    r : unknown
    vartheta_c_in : unknown
    varphi_c_in : unknown
    
    Returns
    -------
    theta_vmec : unknown
    varphi_vmec : unknown
    
    """
    theta_vmec, varphi_vmec = _pysimple.f90wrap_can_to_vmec(r=r, \
        vartheta_c_in=vartheta_c_in, varphi_c_in=varphi_c_in)
    return theta_vmec, varphi_vmec

def deallocate_can_coord():
    """
    deallocate_can_coord()
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 930-938
    
    
    """
    _pysimple.f90wrap_deallocate_can_coord()

def vmec_to_can(r, theta, varphi):
    """
    vartheta_c, varphi_c = vmec_to_can(r, theta, varphi)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/get_canonical_coordinates.f90.i \
        lines 943-1097
    
    Parameters
    ----------
    r : unknown
    theta : unknown
    varphi : unknown
    
    Returns
    -------
    vartheta_c : unknown
    varphi_c : unknown
    
    ------------------------------------------
    """
    vartheta_c, varphi_c = _pysimple.f90wrap_vmec_to_can(r=r, theta=theta, \
        varphi=varphi)
    return vartheta_c, varphi_c

def spline_vmec_data():
    """
    spline_vmec_data()
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 6-298
    
    
    ------------------------------------
     Begin poloidal flux($A_\varphi$):
    """
    _pysimple.f90wrap_spline_vmec_data()

def deallocate_vmec_spline(mode):
    """
    deallocate_vmec_spline(mode)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 303-320
    
    Parameters
    ----------
    mode : int
    
    """
    _pysimple.f90wrap_deallocate_vmec_spline(mode=mode)

def splint_vmec_data(s, theta, varphi):
    """
    a_phi, a_theta, da_phi_ds, da_theta_ds, aiota, r, z, alam, dr_ds, dr_dt, dr_dp, \
        dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp = splint_vmec_data(s, theta, \
        varphi)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 326-477
    
    Parameters
    ----------
    s : unknown
    theta : unknown
    varphi : unknown
    
    Returns
    -------
    a_phi : unknown
    a_theta : unknown
    da_phi_ds : unknown
    da_theta_ds : unknown
    aiota : unknown
    r : unknown
    z : unknown
    alam : unknown
    dr_ds : unknown
    dr_dt : unknown
    dr_dp : unknown
    dz_ds : unknown
    dz_dt : unknown
    dz_dp : unknown
    dl_ds : unknown
    dl_dt : unknown
    dl_dp : unknown
    
    ----------------------------
    """
    a_phi, a_theta, da_phi_ds, da_theta_ds, aiota, r, z, alam, dr_ds, dr_dt, dr_dp, \
        dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp = \
        _pysimple.f90wrap_splint_vmec_data(s=s, theta=theta, varphi=varphi)
    return a_phi, a_theta, da_phi_ds, da_theta_ds, aiota, r, z, alam, dr_ds, dr_dt, \
        dr_dp, dz_ds, dz_dt, dz_dp, dl_ds, dl_dt, dl_dp

def vmec_field(s, theta, varphi, a_theta, a_phi, da_theta_ds, da_phi_ds, aiota, \
    sqg, alam, dl_ds, dl_dt, dl_dp, bctrvr_vartheta, bctrvr_varphi, bcovar_r, \
    bcovar_vartheta, bcovar_varphi):
    """
    vmec_field(s, theta, varphi, a_theta, a_phi, da_theta_ds, da_phi_ds, aiota, sqg, \
        alam, dl_ds, dl_dt, dl_dp, bctrvr_vartheta, bctrvr_varphi, bcovar_r, \
        bcovar_vartheta, bcovar_varphi)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 484-527
    
    Parameters
    ----------
    s : unknown
    theta : unknown
    varphi : unknown
    a_theta : unknown
    a_phi : unknown
    da_theta_ds : unknown
    da_phi_ds : unknown
    aiota : unknown
    sqg : unknown
    alam : unknown
    dl_ds : unknown
    dl_dt : unknown
    dl_dp : unknown
    bctrvr_vartheta : unknown
    bctrvr_varphi : unknown
    bcovar_r : unknown
    bcovar_vartheta : unknown
    bcovar_varphi : unknown
    
    """
    _pysimple.f90wrap_vmec_field(s=s, theta=theta, varphi=varphi, a_theta=a_theta, \
        a_phi=a_phi, da_theta_ds=da_theta_ds, da_phi_ds=da_phi_ds, aiota=aiota, \
        sqg=sqg, alam=alam, dl_ds=dl_ds, dl_dt=dl_dt, dl_dp=dl_dp, \
        bctrvr_vartheta=bctrvr_vartheta, bctrvr_varphi=bctrvr_varphi, \
        bcovar_r=bcovar_r, bcovar_vartheta=bcovar_vartheta, \
        bcovar_varphi=bcovar_varphi)

def splint_iota(s, aiota, daiota_ds):
    """
    splint_iota(s, aiota, daiota_ds)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 532-570
    
    Parameters
    ----------
    s : unknown
    aiota : unknown
    daiota_ds : unknown
    
    """
    _pysimple.f90wrap_splint_iota(s=s, aiota=aiota, daiota_ds=daiota_ds)

def splint_lambda(s, theta, varphi, alam, dl_dt):
    """
    splint_lambda(s, theta, varphi, alam, dl_dt)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 575-646
    
    Parameters
    ----------
    s : unknown
    theta : unknown
    varphi : unknown
    alam : unknown
    dl_dt : unknown
    
    ----------------------------
     Begin interpolation over $\theta$
    """
    _pysimple.f90wrap_splint_lambda(s=s, theta=theta, varphi=varphi, alam=alam, \
        dl_dt=dl_dt)

def s_to_rho_healaxis(m, ns, nrho, nheal, arr_in, arr_out):
    """
    s_to_rho_healaxis(m, ns, nrho, nheal, arr_in, arr_out)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 678-732
    
    Parameters
    ----------
    m : int
    ns : int
    nrho : int
    nheal : int
    arr_in : unknown array
    arr_out : unknown array
    
    """
    _pysimple.f90wrap_s_to_rho_healaxis(m=m, ns=ns, nrho=nrho, nheal=nheal, \
        arr_in=arr_in, arr_out=arr_out)

def determine_nheal_for_axis(m, ns, arr_in):
    """
    nheal = determine_nheal_for_axis(m, ns, arr_in)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 736-800
    
    Parameters
    ----------
    m : int
    ns : int
    arr_in : unknown array
    
    Returns
    -------
    nheal : int
    
    """
    nheal = _pysimple.f90wrap_determine_nheal_for_axis(m=m, ns=ns, arr_in=arr_in)
    return nheal

def volume_and_b00(volume, b00):
    """
    volume_and_b00(volume, b00)
    
    
    Defined at \
        /home/dan/Documents/University/Master/Thesis/repos/SIMPLE/build/spline_vmec_data.f90.i \
        lines 804-852
    
    Parameters
    ----------
    volume : unknown
    b00 : unknown
    
    """
    _pysimple.f90wrap_volume_and_b00(volume=volume, b00=b00)

