#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _CaDynamics_E2_reg(void);
extern void _Ca_HVA_reg(void);
extern void _Ca_LVAst_reg(void);
extern void _epsp_reg(void);
extern void _Ih_reg(void);
extern void _Im_reg(void);
extern void _int2pyr_reg(void);
extern void _K_Pst_reg(void);
extern void _K_Tst_reg(void);
extern void _Nap_Et2_reg(void);
extern void _NaTa_t_reg(void);
extern void _NaTg_reg(void);
extern void _NaTs2_t_reg(void);
extern void _netgaba_reg(void);
extern void _netglutamate_reg(void);
extern void _pyr2pyr_reg(void);
extern void _SK_E2_reg(void);
extern void _SKv3_1_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," modfiles/CaDynamics_E2.mod");
    fprintf(stderr," modfiles/Ca_HVA.mod");
    fprintf(stderr," modfiles/Ca_LVAst.mod");
    fprintf(stderr," modfiles/epsp.mod");
    fprintf(stderr," modfiles/Ih.mod");
    fprintf(stderr," modfiles/Im.mod");
    fprintf(stderr," modfiles/int2pyr.mod");
    fprintf(stderr," modfiles/K_Pst.mod");
    fprintf(stderr," modfiles/K_Tst.mod");
    fprintf(stderr," modfiles/Nap_Et2.mod");
    fprintf(stderr," modfiles/NaTa_t.mod");
    fprintf(stderr," modfiles/NaTg.mod");
    fprintf(stderr," modfiles/NaTs2_t.mod");
    fprintf(stderr," modfiles/netgaba.mod");
    fprintf(stderr," modfiles/netglutamate.mod");
    fprintf(stderr," modfiles/pyr2pyr.mod");
    fprintf(stderr," modfiles/SK_E2.mod");
    fprintf(stderr," modfiles/SKv3_1.mod");
    fprintf(stderr," modfiles/vecevent.mod");
    fprintf(stderr, "\n");
  }
  _CaDynamics_E2_reg();
  _Ca_HVA_reg();
  _Ca_LVAst_reg();
  _epsp_reg();
  _Ih_reg();
  _Im_reg();
  _int2pyr_reg();
  _K_Pst_reg();
  _K_Tst_reg();
  _Nap_Et2_reg();
  _NaTa_t_reg();
  _NaTg_reg();
  _NaTs2_t_reg();
  _netgaba_reg();
  _netglutamate_reg();
  _pyr2pyr_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
  _vecevent_reg();
}
