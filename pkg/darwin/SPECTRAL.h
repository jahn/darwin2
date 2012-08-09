C $Header$
C $Name$

C SPECTRAL.h 
C description: spectral runtime parameters and fields

C  darwin_waves         :: 'central' wavelengths of wavebands (nm)
C  darwin_wavebands     :: waveband boundaries (nm)
      COMMON/SPECTRAL_PARAMS_R/darwin_waves, darwin_wavebands
      _RL darwin_waves(tlam)
      _RL darwin_wavebands(tlam+1)

C WtouEins    :: W to uEin/s conversion factor
      COMMON/darwin_oasim/ WtouEins
      _RL WtouEins(tlam)

C oasim_ed :: spectral direct downwelling irradiance at surface read from file
C oasim_es :: spectral diffuse downwelling irradiance at surface read from file
#ifdef OASIM
         _RL oasim_ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
         _RL oasim_es(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
         COMMON /SPECTRAL_INPUT/ oasim_ed, oasim_es
#endif

      COMMON /SPECTRAL_OUPUT/
#ifdef DAR_DIAG_IRR
     &      Edave,Esave,Euave,Eutave,
#endif
#ifdef DAR_DIAG_IRR_AMPS
     &      c1ave,c2ave,
#endif
#ifdef DAR_DIAG_ABSORP
     &      aave,
#endif
#ifdef DAR_DIAG_SCATTER
     &      btave,
     &      bbave,
#endif
#ifdef DAR_DIAG_PART_SCATTER
     &      apartave,
     &      btpartave,
     &      bbpartave,
#endif
     &      rmudave

#ifdef DAR_DIAG_IRR
       _RL  Edave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  Esave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  Euave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  Eutave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
#endif
#ifdef DAR_DIAG_IRR_AMPS
       _RL  c1ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  c2ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
#endif
#ifdef DAR_DIAG_ABSORP
       _RL  aave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
#endif
#ifdef DAR_DIAG_SCATTER
       _RL  btave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  bbave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
#endif
#ifdef DAR_DIAG_PART_SCATTER
       _RL  apartave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  bbpartave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
       _RL  btpartave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy,tlam)
#endif
       _RL  rmudave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

