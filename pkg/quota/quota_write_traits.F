C $Name:  $

#include "CPP_OPTIONS.h"
#include "DARWIN_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#ifdef ALLOW_DARWIN
#ifdef ALLOW_QUOTA

c ==========================================================
c SUBROUTINE QUOTA_WRITE_TRAITS()
c ==========================================================
c
        SUBROUTINE QUOTA_WRITE_TRAITS(myThid)

        IMPLICIT NONE

#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "QUOTA_SIZE.h"
#include "QUOTA.h"
c#include "GCHEM.h"
#include "DARWIN_IO.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
C     === Local variables ===
      INTEGER IniUnit1

      _RL vmaxDIC(npmax)
      _RL vmaxNO3(npmax)
      _RL vmaxNO2(npmax)
      _RL vmaxNH4(npmax)
      _RL vmaxFeT(npmax)
      _RL ksatNO3(npmax)
      _RL ksatNO2(npmax)
      _RL ksatNH4(npmax)
      _RL ksatFeT(npmax)
      _RL Qnmin(npmax)
      _RL Qfemin(npmax)
      _RL Qnmax(npmax)
      _RL Qfemax(npmax)
      _RL kexcc(npmax)
      _RL kexcn(npmax)
      _RL kexcfe(npmax)
      _RL beta_graz_c(npmax)
      _RL beta_graz_n(npmax)
      _RL beta_graz_fe(npmax)
      _RL beta_mort_c(npmax)
      _RL beta_mort_n(npmax)
      _RL beta_mort_fe(npmax)

      _RL orgsinkD, orgsinkP
      _RL KDOC, KDON, KDOFe
      _RL KPOC, KPON, KPOFe

      NAMELIST/DARWIN_TRAITS/
     &  vmaxDIC,
     &  vmaxNO3,
     &  vmaxNO2,
     &  vmaxNH4,
     &  vmaxFeT,
     &  ksatNO3,
     &  ksatNO2,
     &  ksatNH4,
     &  ksatFeT,
     &  Qnmin,
     &  Qfemin,
     &  Qnmax,
     &  Qfemax,
     &  kexcc,
     &  kexcn,
     &  kexcfe,
     &  beta_graz_c,
     &  beta_graz_n,
     &  beta_graz_fe,
     &  beta_mort_c,
     &  beta_mort_n,
     &  beta_mort_fe,
     &  qcarbon,
     &  alphachl,
     &  pp_opt,
     &  pp_sig,
     &  graz,
     &  kg,
     &  kmort,
     &  respiration,
     &  biosink,
     &  motility,
     &  autotrophy,
     &  phytoTempCoeff,
     &  phytoTempExp1,
     &  phytoTempExp2,
     &  phytoTempOptimum,
     &  phytoDecayPower,
     &  pft,
     &  nsource,
     &  use_NO3,
     &  use_Si,
     &  biovol,
     &  graz_pref,
     &  assim_graz

      NAMELIST/DARWIN_PARAMS/
     &  biosynth,
     &  amminhib,
     &  chl2nmax,
     &  ass_eff,hill,Lambda,
     &  orgsinkD, orgsinkP,
     &  KDOC, KDON, KDOFe,
     &  KPOC, KPON, KPOFe,
     &  amm2nrite,
     &  nrite2nrate,
     &  PAR0,
     &  k_w,
     &  k_chl,
     &  tempcoeff1, tempcoeff2_big,
     &  tempcoeff2_small, tempcoeff3,
     &  tempmax, temprange, tempnorm,tempdecay,
     &  ngrowfac,ilight,
#ifdef IRON_SED_SOURCE
     &  depthfesed, fesedflux, fesedflux_pcm,
#endif
#ifdef PART_SCAV
     &  scav_rat, scav_inter, scav_exp,
#endif
     &  hollingtype,
     &  ns,
     &  ngroups,
     &  k0, kc, parfrac, parconv,
     &  chlpmax, chlpmin, istar

      vmaxDIC(:) = vmaxi(iDIC, :)
      vmaxNO3(:) = vmaxi(iNO3, :)
      vmaxNO2(:) = vmaxi(iNO2, :)
      vmaxNH4(:) = vmaxi(iNH4, :)
      vmaxFeT(:) = vmaxi(iFeT, :)
      ksatNO3(:) = kn(iNO3, :)
      ksatNO2(:) = kn(iNO2, :)
      ksatNH4(:) = kn(iNH4, :)
      ksatFeT(:) = kn(iFeT, :)
      Qnmin(:) = qmin(iNitr, :)
      Qfemin(:) = qmin(iIron, :) 
      Qnmax(:) = qmax(iNitr, :)
      Qfemax(:) = qmax(iIron, :)
      kexcc(:) = kexc(iCarb, :)
      kexcn(:) = kexc(iNitr, :)
      kexcfe(:) = kexc(iIron, :)
      beta_graz_c(:) = beta_graz(iCarb, :)
      beta_graz_n(:) = beta_graz(iNitr, :)
      beta_graz_fe(:) = beta_graz(iIron, :)
      beta_mort_c(:) = beta_mort(iCarb, :)
      beta_mort_n(:) = beta_mort(iNitr, :)
      beta_mort_fe(:) = beta_mort(iIron, :)

      IF ( myProcId.EQ.0 .AND. myThid.EQ.1 ) THEN
        CALL MDSFINDUNIT( IniUnit1, mythid )
        OPEN(IniUnit1,file='darwin_traits.nml',status='unknown')
        WRITE(IniUnit1, nml=darwin_traits)
        CLOSE(IniUnit1)
      ENDIF

      orgsinkD = orgsink(1)
      orgsinkP = orgsink(2)
      KDOC = remin(iCarb, 1)
      KDON = remin(iNitr, 1)
      KDOFe= remin(iIron, 1)
      KPOC = remin(iCarb, 2)
      KPON = remin(iNitr, 2)
      KPOFe= remin(iIron, 2)

      IF ( myProcId.EQ.0 .AND. myThid.EQ.1 ) THEN
        CALL MDSFINDUNIT( IniUnit1, mythid )
        OPEN(IniUnit1,file='darwin_params.nml',status='unknown')
        WRITE(IniUnit1, nml=darwin_params)
        CLOSE(IniUnit1)
      ENDIF

      RETURN
      END
#endif  /*ALLOW_QUOTA*/
#endif  /*ALLOW_DARWIN*/
#endif  /*ALLOW_PTRACERS*/
c ==========================================================
