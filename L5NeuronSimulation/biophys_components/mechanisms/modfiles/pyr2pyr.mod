:Pyramidal Cells to Pyramidal Cells AMPA+NMDA with local Ca2+ pool

NEURON {
	POINT_PROCESS pyr2pyr
	:USEION ca READ eca	
	NONSPECIFIC_CURRENT inmda, iampa
        RANGE tau_r_AMPA, tau_d_AMPA, tau_r_NMDA, tau_d_NMDA
	RANGE initW
	RANGE Cdur_nmda, AlphaTmax_nmda, Beta_nmda, Erev_nmda, gbar_nmda, W_nmda, on_nmda, g_nmda
	RANGE Cdur_ampa, AlphaTmax_ampa, Beta_ampa, Erev_ampa, gbar_ampa, W, on_ampa, g_ampa
	RANGE thr_rp
	RANGE F, f, tauF, D1, d1, tauD1, D2, d2, tauD2
	RANGE facfactor
	RANGE i, iampa, inmda, g_AMPA, g_NMDA, e, initW
	:Release probability
	RANGE random, P, P_0

	THREADSAFE
	POINTER randObjPtr
}

UNITS { 
	(mV) = (millivolt)
        (nA) = (nanoamp)
	(uS) = (microsiemens)
	FARADAY = 96485 (coul)
	pi = 3.141592 (1)
}

PARAMETER {

	:srcid = -1 (1)
	:destid = -1 (1)
	:type = -1
	mggate
	
        tau_r_AMPA = 0.2   (ms)  : dual-exponential conductance profile
        tau_d_AMPA = 1.7    (ms)  : IMPORTANT: tau_r < tau_d
	tau_r_NMDA = 0.29   (ms) : dual-exponential conductance profile
        tau_d_NMDA = 43     (ms) : IMPORTANT: tau_r < tau_d
	
	Cdur_nmda = 16.7650 (ms)
	AlphaTmax_nmda = .2659 (/ms)
	Beta_nmda = 0.008 (/ms)
	Erev_nmda = 0 (mV)
	gbar_nmda = .5e-3 (uS)

	Cdur_ampa = 1.4210 (ms)
	AlphaTmax_ampa = 3.8142 (/ms)
	Beta_ampa =  0.1429(/ms) :0.1429 as original 0.2858 as half,0.07145 as twice
	Erev_ampa = 0 (mV)
	gbar_ampa = 1e-3 (uS)

	mg = 1   (mM)  : initial concentration of mg2+
	:eca = 120

	:Cainf = 50e-6 (mM)
	:pooldiam =  1.8172 (micrometer)
	:z = 2
	:neuroM = 0
	:tauCa = 50 (ms)
	:P0 = .015
	:fCa = .024
	
	:lambda1 = 40 : 60 : 12 :80: 20 : 15 :8 :5: 2.5
	:lambda2 = .03
	:threshold1 = 0.4 :  0.45 : 0.35 :0.35:0.2 :0.50 (uM)
	:threshold2 = 0.55 : 0.50 : 0.40 :0.4 :0.3 :0.60 (uM)

	initW = 5.0 : 1.0 :  0.9 : 0.8 : 2 : 10 : 6 :1.5
	:fmax = 3 : 2.5 : 4 : 2 : 3 : 1.5 : 3
	:fmin = .8
	
	:DAstart1 = 39500
	:DAstop1 = 40000	
	:DAstart2 = 35900
	:DAstop2 = 36000	

	:DA_t1 = 1.2
	:DA_t2 = 0.8 : 0.9
    :DA_t3 = 0.9
	:DA_S = 1.3 : 0.95 : 0.6	
	:Beta1 = 0.001  (/ms) : 1/decay time for neuromodulators
	:Beta2 = 0.0001  (/ms)

	thr_rp = 1 : .7
	
	facfactor = 1
	: the (1) is needed for the range limits to be effective
        f = 0 (1) < 0, 1e9 >    : facilitation
        tauF = 20 (ms) < 1e-9, 1e9 >
        d1 = 0.95 (1) < 0, 1 >     : fast depression
        tauD1 = 40 (ms) < 1e-9, 1e9 >
        d2 = 0.9 (1) < 0, 1 >     : slow depression
        tauD2 = 70 (ms) < 1e-9, 1e9 >	

	P_0 = 1 (1) < 0, 1 >               : base release probability	
}

ASSIGNED {
	v (mV)

	inmda (nA)
	g_NMDA (uS)
	on_nmda
	W_nmda

	iampa (nA)
	g_AMPA (uS)
	on_ampa
	limitW
	i (nA)
	t0 (ms)

	factor_AMPA
	factor_NMDA
	rp
	tsyn
	
	fa
	F
	D1
	D2

	:Release probability
	P				        : instantaneous release probability
    randObjPtr              : pointer to a hoc random number generator Random.uniform(0,1)
    random   
}

STATE {

        A_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_r_AMPA
        B_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_d_AMPA
	A_NMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_r_NMDA
        B_NMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_d_NMDA
}


INITIAL {
        LOCAL tp_AMPA, tp_NMDA
        
	A_AMPA = 0
        B_AMPA = 0
	
	A_NMDA = 0
	B_NMDA = 0
        
	tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA) :time to peak of the conductance
	tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA) :time to peak of the conductance
        
	factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) :AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
        factor_AMPA = 1/factor_AMPA
	
	factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA) :NMDA Normalization factor - so that when t = tp_NMDA, gsyn = gpeak
        factor_NMDA = 1/factor_NMDA
	
	t0 = -1

	fa =0
	F = 1
	D1 = 1
	D2 = 1

	P = P_0
	random = 1
}

BREAKPOINT {

        SOLVE state METHOD cnexp
	mggate = 1 / (1 + exp(0.08  (/mV) * -(v)) * (mg / 3.57 (mM))) :mggate kinetics - Jahr & Stevens 1990
        g_AMPA = facfactor*initW*.001*(B_AMPA-A_AMPA) :compute time varying conductance as the difference of state variables B_AMPA and A_AMPA
	g_NMDA = facfactor*initW*.001*(B_NMDA-A_NMDA) * mggate :compute time varying conductance as the difference of state variables B_NMDA and A_NMDA and mggate kinetics
        iampa = g_AMPA*(v-Erev_ampa) :compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal
	inmda = g_NMDA*(v-Erev_nmda) :compute the NMDA driving force based on the time varying conductance, membrane potential, and NMDA reversal
	i = iampa + inmda
	
}

DERIVATIVE state{

        A_AMPA' = -A_AMPA/tau_r_AMPA
        B_AMPA' = -B_AMPA/tau_d_AMPA
	A_NMDA' = -A_NMDA/tau_r_NMDA
        B_NMDA' = -B_NMDA/tau_d_NMDA
}


NET_RECEIVE (weight,weight_AMPA, weight_NMDA, Pv, Pr, u, tsyn (ms)){
	weight_AMPA = weight
	weight_NMDA = weight
        A_AMPA = A_AMPA + weight_AMPA*factor_AMPA
        B_AMPA = B_AMPA + weight_AMPA*factor_AMPA
	A_NMDA = A_NMDA + weight_NMDA*factor_NMDA
        B_NMDA = B_NMDA + weight_NMDA*factor_NMDA
	
	random = randGen()
	         
	if (flag == 0 && random < P_0) {   : Short term plasticity was implemented(Varela et. al 1997):
		F  = 1 + (F-1)* exp(-(t - tsyn)/tauF)
		D1 = 1 - (1-D1)*exp(-(t - tsyn)/tauD1)
		D2 = 1 - (1-D2)*exp(-(t - tsyn)/tauD2)
		tsyn = t
		
		facfactor = F * D1 * D2	
		F = F * f
		if (F > 30) { 
			F=30	
		}
		if (facfactor < 0.2) { 
			facfactor=0.2
		}	
		D1 = D1 * d1
		D2 = D2 * d2
	}
}

:::::::::::: FUNCTIONs and PROCEDUREs ::::::::::::

FUNCTION sfunc (v (mV)) {
	sfunc = 1 / (1 + exp(0.08  (/mV) * -(v)) * (mg / 3.57 (mM))) :mggate kinetics - Jahr & Stevens 1990
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION randGen() {
VERBATIM
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
ENDVERBATIM
}

PROCEDURE setRandObjRef() {
VERBATIM
   void** pv4 = (void**)(&_p_randObjPtr);
   if (ifarg(1)) {
      *pv4 = nrn_random_arg(1);
   }else{
      *pv4 = (void*)0;
   }
ENDVERBATIM
}
:FUNCTION unirand() {    : uniform random numbers between 0 and 1
:        unirand = scop_random()
:}
