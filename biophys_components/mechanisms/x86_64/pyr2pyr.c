/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__pyr2pyr
#define _nrn_initial _nrn_initial__pyr2pyr
#define nrn_cur _nrn_cur__pyr2pyr
#define _nrn_current _nrn_current__pyr2pyr
#define nrn_jacob _nrn_jacob__pyr2pyr
#define nrn_state _nrn_state__pyr2pyr
#define _net_receive _net_receive__pyr2pyr 
#define release release__pyr2pyr 
#define setRandObjRef setRandObjRef__pyr2pyr 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define initW _p[0]
#define Cdur_nmda _p[1]
#define AlphaTmax_nmda _p[2]
#define Beta_nmda _p[3]
#define Erev_nmda _p[4]
#define gbar_nmda _p[5]
#define Cdur_ampa _p[6]
#define AlphaTmax_ampa _p[7]
#define Beta_ampa _p[8]
#define Erev_ampa _p[9]
#define gbar_ampa _p[10]
#define ECa _p[11]
#define Cainf _p[12]
#define pooldiam _p[13]
#define z _p[14]
#define tauCa _p[15]
#define P0 _p[16]
#define fCa _p[17]
#define lambda1 _p[18]
#define lambda2 _p[19]
#define threshold1 _p[20]
#define threshold2 _p[21]
#define fmax _p[22]
#define fmin _p[23]
#define facfactor _p[24]
#define f _p[25]
#define tauF _p[26]
#define d1 _p[27]
#define tauD1 _p[28]
#define d2 _p[29]
#define tauD2 _p[30]
#define aACH _p[31]
#define bACH _p[32]
#define wACH _p[33]
#define aDA _p[34]
#define bDA _p[35]
#define wDA _p[36]
#define P_0 _p[37]
#define i_nmda _p[38]
#define g_nmda _p[39]
#define on_nmda _p[40]
#define W_nmda _p[41]
#define i_ampa _p[42]
#define g_ampa _p[43]
#define on_ampa _p[44]
#define W_ampa _p[45]
#define ICa _p[46]
#define iCatotal _p[47]
#define Wmax _p[48]
#define Wmin _p[49]
#define maxChange _p[50]
#define normW _p[51]
#define scaleW _p[52]
#define pregid _p[53]
#define postgid _p[54]
#define calcium _p[55]
#define F _p[56]
#define D1 _p[57]
#define D2 _p[58]
#define P _p[59]
#define random _p[60]
#define r_nmda _p[61]
#define r_ampa _p[62]
#define Capoolcon _p[63]
#define t0 _p[64]
#define Afactor _p[65]
#define dW_ampa _p[66]
#define tsyn _p[67]
#define fa _p[68]
#define Dr_nmda _p[69]
#define Dr_ampa _p[70]
#define DCapoolcon _p[71]
#define v _p[72]
#define _g _p[73]
#define _tsav _p[74]
#define _nd_area  *_ppvar[0]._pval
#define randObjPtr	*_ppvar[2]._pval
#define _p_randObjPtr	_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  2;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_eta();
 static double _hoc_omega();
 static double _hoc_randGen();
 static double _hoc_setRandObjRef();
 static double _hoc_sfunc();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "eta", _hoc_eta,
 "omega", _hoc_omega,
 "randGen", _hoc_randGen,
 "setRandObjRef", _hoc_setRandObjRef,
 "sfunc", _hoc_sfunc,
 0, 0
};
#define eta eta_pyr2pyr
#define omega omega_pyr2pyr
#define randGen randGen_pyr2pyr
#define sfunc sfunc_pyr2pyr
 extern double eta( _threadargsprotocomma_ double );
 extern double omega( _threadargsprotocomma_ double , double , double );
 extern double randGen( _threadargsproto_ );
 extern double sfunc( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define ACH ACH_pyr2pyr
 double ACH = 1;
#define DA DA_pyr2pyr
 double DA = 1;
#define LearningShutDown LearningShutDown_pyr2pyr
 double LearningShutDown = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "P_0", 0, 1,
 "d2", 0, 1,
 "d1", 0, 1,
 "f", 0, 1e+09,
 "tauD2", 1e-09, 1e+09,
 "tauD1", 1e-09, 1e+09,
 "tauF", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Cdur_nmda", "ms",
 "AlphaTmax_nmda", "/ms",
 "Beta_nmda", "/ms",
 "Erev_nmda", "mV",
 "gbar_nmda", "uS",
 "Cdur_ampa", "ms",
 "AlphaTmax_ampa", "/ms",
 "Beta_ampa", "/ms",
 "Erev_ampa", "mV",
 "gbar_ampa", "uS",
 "Cainf", "mM",
 "pooldiam", "micrometer",
 "tauCa", "ms",
 "threshold1", "uM",
 "threshold2", "uM",
 "f", "1",
 "tauF", "ms",
 "d1", "1",
 "tauD1", "ms",
 "d2", "1",
 "tauD2", "ms",
 "P_0", "1",
 "i_nmda", "nA",
 "g_nmda", "uS",
 "i_ampa", "nA",
 "g_ampa", "uS",
 "ICa", "mA",
 "iCatotal", "mA",
 0,0
};
 static double Capoolcon0 = 0;
 static double delta_t = 0.01;
 static double r_ampa0 = 0;
 static double r_nmda0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ACH_pyr2pyr", &ACH_pyr2pyr,
 "DA_pyr2pyr", &DA_pyr2pyr,
 "LearningShutDown_pyr2pyr", &LearningShutDown_pyr2pyr,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"pyr2pyr",
 "initW",
 "Cdur_nmda",
 "AlphaTmax_nmda",
 "Beta_nmda",
 "Erev_nmda",
 "gbar_nmda",
 "Cdur_ampa",
 "AlphaTmax_ampa",
 "Beta_ampa",
 "Erev_ampa",
 "gbar_ampa",
 "ECa",
 "Cainf",
 "pooldiam",
 "z",
 "tauCa",
 "P0",
 "fCa",
 "lambda1",
 "lambda2",
 "threshold1",
 "threshold2",
 "fmax",
 "fmin",
 "facfactor",
 "f",
 "tauF",
 "d1",
 "tauD1",
 "d2",
 "tauD2",
 "aACH",
 "bACH",
 "wACH",
 "aDA",
 "bDA",
 "wDA",
 "P_0",
 0,
 "i_nmda",
 "g_nmda",
 "on_nmda",
 "W_nmda",
 "i_ampa",
 "g_ampa",
 "on_ampa",
 "W_ampa",
 "ICa",
 "iCatotal",
 "Wmax",
 "Wmin",
 "maxChange",
 "normW",
 "scaleW",
 "pregid",
 "postgid",
 "calcium",
 "F",
 "D1",
 "D2",
 "P",
 "random",
 0,
 "r_nmda",
 "r_ampa",
 "Capoolcon",
 0,
 "randObjPtr",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 75, _prop);
 	/*initialize range parameters*/
 	initW = 5;
 	Cdur_nmda = 17.58;
 	AlphaTmax_nmda = 0.08;
 	Beta_nmda = 0.008;
 	Erev_nmda = 0;
 	gbar_nmda = 0.0006;
 	Cdur_ampa = 5.31;
 	AlphaTmax_ampa = 0.117;
 	Beta_ampa = 0.09;
 	Erev_ampa = 0;
 	gbar_ampa = 0.0017;
 	ECa = 120;
 	Cainf = 5e-05;
 	pooldiam = 1.8172;
 	z = 2;
 	tauCa = 50;
 	P0 = 0.015;
 	fCa = 0.024;
 	lambda1 = 2.5;
 	lambda2 = 0.01;
 	threshold1 = 0.2;
 	threshold2 = 0.4;
 	fmax = 3;
 	fmin = 0.8;
 	facfactor = 1;
 	f = 1;
 	tauF = 1;
 	d1 = 1;
 	tauD1 = 1;
 	d2 = 1;
 	tauD2 = 1;
 	aACH = 1;
 	bACH = 0;
 	wACH = 0;
 	aDA = 1;
 	bDA = 0;
 	wDA = 0;
 	P_0 = 1;
  }
 	_prop->param = _p;
 	_prop->param_size = 75;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _pyr2pyr_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 75, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 pyr2pyr /home/mizzou/Internship/L5Neuron/GoodL5NeuronSimulation/biophys_components/mechanisms/x86_64/pyr2pyr.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.0;
 static double pi = 3.141592;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int setRandObjRef(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int release(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   if ( t0 > 0.0 ) {
     if ( t - t0 < Cdur_nmda ) {
       on_nmda = 1.0 ;
       }
     else {
       on_nmda = 0.0 ;
       }
     if ( t - t0 < Cdur_ampa ) {
       on_ampa = 1.0 ;
       }
     else {
       on_ampa = 0.0 ;
       }
     }
   Dr_nmda = AlphaTmax_nmda * on_nmda * ( 1.0 - r_nmda ) - Beta_nmda * r_nmda ;
   Dr_ampa = AlphaTmax_ampa * on_ampa * ( 1.0 - r_ampa ) - Beta_ampa * r_ampa ;
   dW_ampa = eta ( _threadargscomma_ Capoolcon ) * ( lambda1 * omega ( _threadargscomma_ Capoolcon , threshold1 , threshold2 ) - lambda2 * W_ampa ) * dt ;
   if ( fabs ( dW_ampa ) > maxChange ) {
     if ( dW_ampa < 0.0 ) {
       dW_ampa = - 1.0 * maxChange ;
       }
     else {
       dW_ampa = maxChange ;
       }
     }
   normW = ( W_ampa - Wmin ) / ( Wmax - Wmin ) ;
   if ( dW_ampa < 0.0 ) {
     scaleW = sqrt ( fabs ( normW ) ) ;
     }
   else {
     scaleW = sqrt ( fabs ( 1.0 - normW ) ) ;
     }
   W_ampa = W_ampa + dW_ampa * scaleW * ( 1.0 + ( wACH * ( ACH - 1.0 ) ) ) * LearningShutDown ;
   if ( W_ampa > Wmax ) {
     W_ampa = Wmax ;
     }
   else if ( W_ampa < Wmin ) {
     W_ampa = Wmin ;
     }
   g_nmda = gbar_nmda * r_nmda * facfactor ;
   i_nmda = W_nmda * g_nmda * ( v - Erev_nmda ) * sfunc ( _threadargscomma_ v ) ;
   g_ampa = gbar_ampa * r_ampa * facfactor ;
   i_ampa = W_ampa * g_ampa * ( v - Erev_ampa ) * ( 1.0 + ( bACH * ( ACH - 1.0 ) ) ) * ( aDA + ( bDA * ( DA - 1.0 ) ) ) ;
   ICa = P0 * g_nmda * ( v - ECa ) * sfunc ( _threadargscomma_ v ) ;
   DCapoolcon = - fCa * Afactor * ICa + ( Cainf - Capoolcon ) / tauCa ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 if ( t0 > 0.0 ) {
   if ( t - t0 < Cdur_nmda ) {
     on_nmda = 1.0 ;
     }
   else {
     on_nmda = 0.0 ;
     }
   if ( t - t0 < Cdur_ampa ) {
     on_ampa = 1.0 ;
     }
   else {
     on_ampa = 0.0 ;
     }
   }
 Dr_nmda = Dr_nmda  / (1. - dt*( ( AlphaTmax_nmda * on_nmda )*( ( ( - 1.0 ) ) ) - ( Beta_nmda )*( 1.0 ) )) ;
 Dr_ampa = Dr_ampa  / (1. - dt*( ( AlphaTmax_ampa * on_ampa )*( ( ( - 1.0 ) ) ) - ( Beta_ampa )*( 1.0 ) )) ;
 dW_ampa = eta ( _threadargscomma_ Capoolcon ) * ( lambda1 * omega ( _threadargscomma_ Capoolcon , threshold1 , threshold2 ) - lambda2 * W_ampa ) * dt ;
 if ( fabs ( dW_ampa ) > maxChange ) {
   if ( dW_ampa < 0.0 ) {
     dW_ampa = - 1.0 * maxChange ;
     }
   else {
     dW_ampa = maxChange ;
     }
   }
 normW = ( W_ampa - Wmin ) / ( Wmax - Wmin ) ;
 if ( dW_ampa < 0.0 ) {
   scaleW = sqrt ( fabs ( normW ) ) ;
   }
 else {
   scaleW = sqrt ( fabs ( 1.0 - normW ) ) ;
   }
 W_ampa = W_ampa + dW_ampa * scaleW * ( 1.0 + ( wACH * ( ACH - 1.0 ) ) ) * LearningShutDown ;
 if ( W_ampa > Wmax ) {
   W_ampa = Wmax ;
   }
 else if ( W_ampa < Wmin ) {
   W_ampa = Wmin ;
   }
 g_nmda = gbar_nmda * r_nmda * facfactor ;
 i_nmda = W_nmda * g_nmda * ( v - Erev_nmda ) * sfunc ( _threadargscomma_ v ) ;
 g_ampa = gbar_ampa * r_ampa * facfactor ;
 i_ampa = W_ampa * g_ampa * ( v - Erev_ampa ) * ( 1.0 + ( bACH * ( ACH - 1.0 ) ) ) * ( aDA + ( bDA * ( DA - 1.0 ) ) ) ;
 ICa = P0 * g_nmda * ( v - ECa ) * sfunc ( _threadargscomma_ v ) ;
 DCapoolcon = DCapoolcon  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauCa )) ;
  return 0;
}
 /*END CVODE*/
 static int release (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   if ( t0 > 0.0 ) {
     if ( t - t0 < Cdur_nmda ) {
       on_nmda = 1.0 ;
       }
     else {
       on_nmda = 0.0 ;
       }
     if ( t - t0 < Cdur_ampa ) {
       on_ampa = 1.0 ;
       }
     else {
       on_ampa = 0.0 ;
       }
     }
    r_nmda = r_nmda + (1. - exp(dt*(( AlphaTmax_nmda * on_nmda )*( ( ( - 1.0 ) ) ) - ( Beta_nmda )*( 1.0 ))))*(- ( ( ( AlphaTmax_nmda )*( on_nmda ) )*( ( 1.0 ) ) ) / ( ( ( AlphaTmax_nmda )*( on_nmda ) )*( ( ( - 1.0 ) ) ) - ( Beta_nmda )*( 1.0 ) ) - r_nmda) ;
    r_ampa = r_ampa + (1. - exp(dt*(( AlphaTmax_ampa * on_ampa )*( ( ( - 1.0 ) ) ) - ( Beta_ampa )*( 1.0 ))))*(- ( ( ( AlphaTmax_ampa )*( on_ampa ) )*( ( 1.0 ) ) ) / ( ( ( AlphaTmax_ampa )*( on_ampa ) )*( ( ( - 1.0 ) ) ) - ( Beta_ampa )*( 1.0 ) ) - r_ampa) ;
   dW_ampa = eta ( _threadargscomma_ Capoolcon ) * ( lambda1 * omega ( _threadargscomma_ Capoolcon , threshold1 , threshold2 ) - lambda2 * W_ampa ) * dt ;
   if ( fabs ( dW_ampa ) > maxChange ) {
     if ( dW_ampa < 0.0 ) {
       dW_ampa = - 1.0 * maxChange ;
       }
     else {
       dW_ampa = maxChange ;
       }
     }
   normW = ( W_ampa - Wmin ) / ( Wmax - Wmin ) ;
   if ( dW_ampa < 0.0 ) {
     scaleW = sqrt ( fabs ( normW ) ) ;
     }
   else {
     scaleW = sqrt ( fabs ( 1.0 - normW ) ) ;
     }
   W_ampa = W_ampa + dW_ampa * scaleW * ( 1.0 + ( wACH * ( ACH - 1.0 ) ) ) * LearningShutDown ;
   if ( W_ampa > Wmax ) {
     W_ampa = Wmax ;
     }
   else if ( W_ampa < Wmin ) {
     W_ampa = Wmin ;
     }
   g_nmda = gbar_nmda * r_nmda * facfactor ;
   i_nmda = W_nmda * g_nmda * ( v - Erev_nmda ) * sfunc ( _threadargscomma_ v ) ;
   g_ampa = gbar_ampa * r_ampa * facfactor ;
   i_ampa = W_ampa * g_ampa * ( v - Erev_ampa ) * ( 1.0 + ( bACH * ( ACH - 1.0 ) ) ) * ( aDA + ( bDA * ( DA - 1.0 ) ) ) ;
   ICa = P0 * g_nmda * ( v - ECa ) * sfunc ( _threadargscomma_ v ) ;
    Capoolcon = Capoolcon + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauCa)))*(- ( ( ( - fCa )*( Afactor ) )*( ICa ) + ( ( Cainf ) ) / tauCa ) / ( ( ( ( - 1.0 ) ) ) / tauCa ) - Capoolcon) ;
   }
  return 0;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   random = randGen ( _threadargs_ ) ;
   if ( random < P_0 ) {
     t0 = t ;
     F = 1.0 + ( F - 1.0 ) * exp ( - ( t - tsyn ) / tauF ) ;
     D1 = 1.0 - ( 1.0 - D1 ) * exp ( - ( t - tsyn ) / tauD1 ) ;
     D2 = 1.0 - ( 1.0 - D2 ) * exp ( - ( t - tsyn ) / tauD2 ) ;
     tsyn = t ;
     facfactor = F * D1 * D2 ;
     F = F * f ;
     if ( F > 30.0 ) {
       F = 30.0 ;
       }
     D1 = D1 * d1 ;
     D2 = D2 * d2 ;
     }
   } }
 
double sfunc ( _threadargsprotocomma_ double _lv ) {
   double _lsfunc;
  _lsfunc = 1.0 / ( 1.0 + 0.33 * exp ( - 0.06 * _lv ) ) ;
    
return _lsfunc;
 }
 
static double _hoc_sfunc(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  sfunc ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
double eta ( _threadargsprotocomma_ double _lCani ) {
   double _leta;
 double _ltaulearn , _lP1 , _lP2 , _lP4 , _lCacon ;
 _lP1 = 0.1 ;
   _lP2 = _lP1 * 1e-4 ;
   _lP4 = 1.0 ;
   _lCacon = _lCani * 1e3 ;
   _ltaulearn = _lP1 / ( _lP2 + _lCacon * _lCacon * _lCacon ) + _lP4 ;
   _leta = 1.0 / _ltaulearn * 0.001 ;
   
return _leta;
 }
 
static double _hoc_eta(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  eta ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
double omega ( _threadargsprotocomma_ double _lCani , double _lthreshold1 , double _lthreshold2 ) {
   double _lomega;
 double _lr , _lmid , _lCacon ;
 _lCacon = _lCani * 1e3 ;
   _lr = ( _lthreshold2 - _lthreshold1 ) / 2.0 ;
   _lmid = ( _lthreshold1 + _lthreshold2 ) / 2.0 ;
   if ( _lCacon <= _lthreshold1 ) {
     _lomega = 0.0 ;
     }
   else if ( _lCacon >= _lthreshold2 ) {
     _lomega = 1.0 / ( 1.0 + 50.0 * exp ( - 50.0 * ( _lCacon - _lthreshold2 ) ) ) ;
     }
   else {
     _lomega = - sqrt ( _lr * _lr - ( _lCacon - _lmid ) * ( _lCacon - _lmid ) ) ;
     }
   
return _lomega;
 }
 
static double _hoc_omega(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  omega ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) );
 return(_r);
}
 
/*VERBATIM*/
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
 
double randGen ( _threadargsproto_ ) {
   double _lrandGen;
 
/*VERBATIM*/
   if (_p_randObjPtr) {
      /*
      :Supports separate independent but reproducible streams for
      : each instance. However, the corresponding hoc Random
      : distribution MUST be set to Random.uniform(0,1)
      */
      _lrandGen = nrn_random_pick(_p_randObjPtr);
   }else{
      hoc_execerror("Random object ref not set correctly for randObjPtr"," only via hoc Random");
   }
 
return _lrandGen;
 }
 
static double _hoc_randGen(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  randGen ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
static int  setRandObjRef ( _threadargsproto_ ) {
   
/*VERBATIM*/
   void** pv4 = (void**)(&_p_randObjPtr);
   if (ifarg(1)) {
      *pv4 = nrn_random_arg(1);
   }else{
      *pv4 = (void*)0;
   }
  return 0; }
 
static double _hoc_setRandObjRef(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (_NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRandObjRef ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  Capoolcon = Capoolcon0;
  r_ampa = r_ampa0;
  r_nmda = r_nmda0;
 {
   on_nmda = 0.0 ;
   r_nmda = 0.0 ;
   W_nmda = initW ;
   on_ampa = 0.0 ;
   r_ampa = 0.0 ;
   W_ampa = initW ;
   t0 = - 1.0 ;
   maxChange = ( Wmax - Wmin ) / 10.0 ;
   dW_ampa = 0.0 ;
   Capoolcon = Cainf ;
   Afactor = 1.0 / ( z * FARADAY * 4.0 / 3.0 * pi * pow( ( pooldiam / 2.0 ) , 3.0 ) ) * ( 1e6 ) ;
   tsyn = - 1e30 ;
   fa = 0.0 ;
   F = 1.0 ;
   D1 = 1.0 ;
   D2 = 1.0 ;
   P = P_0 ;
   random = 1.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   }
 _current += i_nmda;
 _current += i_ampa;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   release(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(r_nmda) - _p;  _dlist1[0] = &(Dr_nmda) - _p;
 _slist1[1] = &(r_ampa) - _p;  _dlist1[1] = &(Dr_ampa) - _p;
 _slist1[2] = &(Capoolcon) - _p;  _dlist1[2] = &(DCapoolcon) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/mizzou/Internship/L5Neuron/GoodL5NeuronSimulation/biophys_components/mechanisms/modfiles/pyr2pyr.mod";
static const char* nmodl_file_text = 
  "NEURON {\n"
  "	POINT_PROCESS pyr2pyr\n"
  "	NONSPECIFIC_CURRENT i_nmda, i_ampa\n"
  "	RANGE initW\n"
  "	RANGE Cdur_nmda, AlphaTmax_nmda, Beta_nmda, Erev_nmda, gbar_nmda, W_nmda, on_nmda, g_nmda\n"
  "	RANGE Cdur_ampa, AlphaTmax_ampa, Beta_ampa, Erev_ampa, gbar_ampa, W_ampa, on_ampa, g_ampa\n"
  "	RANGE ECa, ICa, P0, fCa, tauCa, iCatotal\n"
  "	RANGE Cainf, pooldiam, z\n"
  "	RANGE lambda1, lambda2, threshold1, threshold2\n"
  "	RANGE fmax, fmin, Wmax, Wmin, maxChange, normW, scaleW\n"
  "	RANGE pregid,postgid\n"
  "	\n"
  "	:Added by Ali\n"
  "	RANGE F, f, tauF, D1, d1, tauD1, D2, d2, tauD2\n"
  "	RANGE facfactor\n"
  "	RANGE aACH, bACH, aDA, bDA, wACH, wDA, calcium\n"
  "\n"
  "	:Release probability\n"
  "	RANGE random, P, P_0\n"
  "\n"
  "	THREADSAFE\n"
  "	POINTER randObjPtr\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "        (nA) = (nanoamp)\n"
  "	(uS) = (microsiemens)\n"
  "	FARADAY = 96485 (coul)\n"
  "	pi = 3.141592 (1)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  ": parameters are vars assigned by user or changed by hoc. THey appear in nrnpointmenu\n"
  "	initW = 5\n"
  "\n"
  "	Cdur_nmda = 17.58 (ms)\n"
  "	AlphaTmax_nmda = .08 (/ms)\n"
  "	Beta_nmda = 0.008 (/ms)\n"
  "	Erev_nmda = 0 (mV)\n"
  "	gbar_nmda = .6e-3 (uS)\n"
  "\n"
  "	Cdur_ampa = 5.31 (ms)\n"
  "	AlphaTmax_ampa = 0.117 (/ms)\n"
  "	Beta_ampa = 0.090 (/ms)\n"
  "	Erev_ampa = 0 (mV)\n"
  "	gbar_ampa = 1.7e-3 (uS)\n"
  "\n"
  "	ECa = 120\n"
  "\n"
  "	Cainf = 50e-6 (mM)\n"
  "	pooldiam =  1.8172 (micrometer)\n"
  "	z = 2\n"
  "\n"
  "	tauCa = 50 (ms)\n"
  "	P0 = .015\n"
  "	fCa = .024\n"
  "\n"
  "	lambda1 = 2.5\n"
  "	lambda2 = .01\n"
  "	threshold1 = 0.2 (uM)\n"
  "	threshold2 = 0.4 (uM)\n"
  "\n"
  "	fmax = 3\n"
  "	fmin = .8\n"
  "\n"
  "	:Added by Ali\n"
  "	ACH = 1\n"
  "	DA = 1\n"
  "	LearningShutDown = 0\n"
  "\n"
  "	facfactor = 1\n"
  "	: the (1) is needed for the range limits to be effective\n"
  "        f = 1 (1) < 0, 1e9 >    : facilitation\n"
  "        tauF = 1 (ms) < 1e-9, 1e9 >\n"
  "        d1 = 1 (1) < 0, 1 >     : fast depression\n"
  "        tauD1 = 1 (ms) < 1e-9, 1e9 >\n"
  "        d2 = 1 (1) < 0, 1 >     : slow depression\n"
  "        tauD2 = 1 (ms) < 1e-9, 1e9 >\n"
  "		\n"
  "	aACH = 1\n"
  "	bACH = 0\n"
  "	wACH = 0\n"
  "	aDA = 1\n"
  "	bDA = 0\n"
  "	wDA = 0\n"
  "\n"
  "	P_0 = 1 (1) < 0, 1 >               : base release probability\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  ": These are vars calculated by Neuron hoc or by the mechanism mod itself\n"
  "	v (mV)\n"
  "\n"
  "	i_nmda (nA)\n"
  "	g_nmda (uS)\n"
  "	on_nmda\n"
  "	W_nmda\n"
  "\n"
  "	i_ampa (nA)\n"
  "	g_ampa (uS)\n"
  "	on_ampa\n"
  "	W_ampa\n"
  "\n"
  "	t0 (ms)\n"
  "\n"
  "	ICa (mA)\n"
  "	Afactor	(mM/ms/nA)\n"
  "	iCatotal (mA)\n"
  "\n"
  "	dW_ampa\n"
  "	Wmax\n"
  "	Wmin\n"
  "	maxChange\n"
  "	normW\n"
  "	scaleW\n"
  "	\n"
  "	pregid\n"
  "	postgid\n"
  "	\n"
  "	:Added by Ali\n"
  "		calcium\n"
  "\n"
  "		tsyn\n"
  "	\n"
  "		fa\n"
  "		F\n"
  "		D1\n"
  "		D2\n"
  "\n"
  "	:Release probability\n"
  "		P				        : instantaneous release probability\n"
  "    randObjPtr              : pointer to a hoc random number generator Random.uniform(0,1)\n"
  "    random                  : individual instance of random number\n"
  "}\n"
  "\n"
  "STATE { r_nmda r_ampa Capoolcon }\n"
  "\n"
  "INITIAL {\n"
  "	on_nmda = 0\n"
  "	r_nmda = 0\n"
  "	W_nmda = initW\n"
  "\n"
  "	on_ampa = 0\n"
  "	r_ampa = 0\n"
  "	W_ampa = initW\n"
  "\n"
  "	t0 = -1\n"
  "\n"
  "	:Wmax = 2*initW\n"
  "	:Wmin = 0.25*initW\n"
  "	maxChange = (Wmax-Wmin)/10\n"
  "	dW_ampa = 0\n"
  "\n"
  "	Capoolcon = Cainf\n"
  "	Afactor	= 1/(z*FARADAY*4/3*pi*(pooldiam/2)^3)*(1e6)\n"
  "	\n"
  "	:Added by Ali\n"
  "\n"
  "		tsyn = -1e30\n"
  "\n"
  "	fa =0\n"
  "	F = 1\n"
  "	D1 = 1\n"
  "	D2 = 1\n"
  "\n"
  "	P = P_0\n"
  "	random = 1\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE release METHOD cnexp\n"
  "}\n"
  "\n"
  "DERIVATIVE release {\n"
  "	if (t0>0) {\n"
  "		if (t-t0 < Cdur_nmda) {\n"
  "			on_nmda = 1\n"
  "		} else {\n"
  "			on_nmda = 0\n"
  "		}\n"
  "		if (t-t0 < Cdur_ampa) {\n"
  "			on_ampa = 1\n"
  "		} else {\n"
  "			on_ampa = 0\n"
  "		}\n"
  "	}\n"
  "	r_nmda' = AlphaTmax_nmda*on_nmda*(1-r_nmda) -Beta_nmda*r_nmda\n"
  "	r_ampa' = AlphaTmax_ampa*on_ampa*(1-r_ampa) -Beta_ampa*r_ampa\n"
  "\n"
  "	dW_ampa = eta(Capoolcon)*(lambda1*omega(Capoolcon, threshold1, threshold2)-lambda2*W_ampa)*dt\n"
  "\n"
  "	: Limit for extreme large weight changes\n"
  "	if (fabs(dW_ampa) > maxChange) {\n"
  "		if (dW_ampa < 0) {\n"
  "			dW_ampa = -1*maxChange\n"
  "		} else {\n"
  "			dW_ampa = maxChange\n"
  "		}\n"
  "	}\n"
  "\n"
  "	:Normalize the weight change\n"
  "	normW = (W_ampa-Wmin)/(Wmax-Wmin)\n"
  "	if (dW_ampa < 0) {\n"
  "		scaleW = sqrt(fabs(normW))\n"
  "	} else {\n"
  "		scaleW = sqrt(fabs(1.0-normW))\n"
  "	}\n"
  "\n"
  "	W_ampa = W_ampa + dW_ampa*scaleW *(1+ (wACH * (ACH - 1))) * LearningShutDown\n"
  "	\n"
  "	:Weight value limits\n"
  "	if (W_ampa > Wmax) { \n"
  "		W_ampa = Wmax\n"
  "	} else if (W_ampa < Wmin) {\n"
  " 		W_ampa = Wmin\n"
  "	}\n"
  "\n"
  "	g_nmda = gbar_nmda*r_nmda * facfactor\n"
  "	i_nmda = W_nmda*g_nmda*(v - Erev_nmda)*sfunc(v)\n"
  "\n"
  "	g_ampa = gbar_ampa*r_ampa * facfactor\n"
  "	i_ampa = W_ampa*g_ampa*(v - Erev_ampa)  * (1 + (bACH * (ACH-1)))*(aDA + (bDA * (DA-1))) \n"
  "\n"
  "	ICa = P0*g_nmda*(v - ECa)*sfunc(v)\n"
  "	Capoolcon'= -fCa*Afactor*ICa + (Cainf-Capoolcon)/tauCa\n"
  "}\n"
  "NET_RECEIVE(dummy_weight) {\n"
  "	random = randGen()\n"
  "	if (random < P_0)\n"
  "	{\n"
  "		t0 = t :spike time for conductance opening\n"
  "		\n"
  "		:Added by Ali, Synaptic facilitation\n"
  "		F  = 1 + (F-1)* exp(-(t - tsyn)/tauF)\n"
  "		D1 = 1 - (1-D1)*exp(-(t - tsyn)/tauD1)\n"
  "		D2 = 1 - (1-D2)*exp(-(t - tsyn)/tauD2)\n"
  "	:printf(\"%g\\t%g\\t%g\\t%g\\t%g\\t%g\\n\", t, t-tsyn, F, D1, D2, facfactor)\n"
  "		:if (P_0*F*D1 > 1) {\n"
  "		:	P = 1\n"
  "		:} else {\n"
  "		:	P = P_0*F*D1*D2\n"
  "		:}\n"
  "		:random = randGen()\n"
  "		:if (random <= P) {\n"
  "			:net_event(t)\n"
  "		:}\n"
  "\n"
  "		tsyn = t\n"
  "		\n"
  "		facfactor = F * D1 * D2\n"
  "		F = F * f\n"
  "		\n"
  "		if (F > 30) { \n"
  "		F=30\n"
  "		}\n"
  "		D1 = D1 * d1\n"
  "		D2 = D2 * d2\n"
  "	}\n"
  ":printf(\"\\t%g\\t%g\\t%g\\n\", F, D1, D2)\n"
  "	\n"
  "}\n"
  ":::::::::::: FUNCTIONs and PROCEDUREs ::::::::::::\n"
  "FUNCTION sfunc (v (mV)) {\n"
  "	UNITSOFF\n"
  "	sfunc = 1/(1+0.33*exp(-0.06*v))\n"
  "	UNITSON\n"
  "}\n"
  "FUNCTION eta(Cani (mM)) {\n"
  "	LOCAL taulearn, P1, P2, P4, Cacon\n"
  "	P1 = 0.1\n"
  "	P2 = P1*1e-4\n"
  "	P4 = 1\n"
  "	Cacon = Cani*1e3\n"
  "	taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4\n"
  "	eta = 1/taulearn*0.001\n"
  "}\n"
  "FUNCTION omega(Cani (mM), threshold1 (uM), threshold2 (uM)) {\n"
  "	LOCAL r, mid, Cacon\n"
  "	Cacon = Cani*1e3\n"
  "	r = (threshold2-threshold1)/2\n"
  "	mid = (threshold1+threshold2)/2\n"
  "	if (Cacon <= threshold1) { omega = 0}\n"
  "	else if (Cacon >= threshold2) {	omega = 1/(1+50*exp(-50*(Cacon-threshold2)))}\n"
  "	else {omega = -sqrt(r*r-(Cacon-mid)*(Cacon-mid))}\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "double nrn_random_pick(void* r);\n"
  "void* nrn_random_arg(int argpos);\n"
  "ENDVERBATIM\n"
  "\n"
  "FUNCTION randGen() {\n"
  "VERBATIM\n"
  "   if (_p_randObjPtr) {\n"
  "      /*\n"
  "      :Supports separate independent but reproducible streams for\n"
  "      : each instance. However, the corresponding hoc Random\n"
  "      : distribution MUST be set to Random.uniform(0,1)\n"
  "      */\n"
  "      _lrandGen = nrn_random_pick(_p_randObjPtr);\n"
  "   }else{\n"
  "      hoc_execerror(\"Random object ref not set correctly for randObjPtr\",\" only via hoc Random\");\n"
  "   }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE setRandObjRef() {\n"
  "VERBATIM\n"
  "   void** pv4 = (void**)(&_p_randObjPtr);\n"
  "   if (ifarg(1)) {\n"
  "      *pv4 = nrn_random_arg(1);\n"
  "   }else{\n"
  "      *pv4 = (void*)0;\n"
  "   }\n"
  "ENDVERBATIM\n"
  "}\n"
  ;
#endif
