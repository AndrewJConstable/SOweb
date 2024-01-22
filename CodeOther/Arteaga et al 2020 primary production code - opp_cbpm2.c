#include "opp_def.h"
#include "opp_lib.h"

#include <stdio.h>
#include <math.h>

/*--------------------------------------------------------------------------*/

double opp_cbpm2( double chl,
		  double bbp,
		  double irr,
		  double k490,
		  double mld,
		  double zno3,
		  double daylength) {
/*

   !Description:     opp_cbpm2 - computes daily primary productivity using a chl:Carbon ratio.  
                     This is a spectrally resolved version of the cbpm, using nine separate
                     wavelengths.  It is also depth resolved, integrating the effects from
                     the surface down to a fixed depth of 200 m.

                     The cbpm2 algorithm estimates productivity using chl (m-1), bbp (m-1), 
		     surface irradiance (Einsteins m-2 d-1), k490 (m-1), mld (m), zno3 (m) 
                     and day length (hours).  

Net primary productivity is carbon * growth rate, where carbon is proportional to particulate
backscatter

	carbon = 13000 * (bbp - 0.00035)

and growth rate is a function of nutrient and temperature stress (f(nut,T) and photoacclimation 
(f(Ig))

	growth rate (u) = umax * f(nut,T) * f(Ig)

where:

	umax = 2

	f(nut,T) = ((Chl/C)sat - y0) / ((Chl/C)max - y0)

	f(Ig) = 1 - exp (-5 * Ig)


and: 

	(Chl/C)sat = ratio of satellite observed chl and carbon (carbon from bbp)

	(Chl/C)max = 0.022 + (0.045-0.022) * exp (-3 * Ig)

	Ig = median mixed layer light level 
	   = surface irradiance * exp (-k(lambda) * MLD/2)

The above items are analyzed for nine separate wavelengths, and is vertically resolved to a depth
of 200 m.

For more details, please see the paper by Westberry, et al (2008)


   !Input Parameters:  
      chl            chlorophyll concentration
      bbp            backscatter
      irr            Photosynthetically available radiation in Einsteins per
                     day per square meter
      k490           absorbence at 490nm
      mld            mixing layer depth in meters
      zno3           depth of the nitrocline
      daylength      length of the day in decimal hours.

   !Output Parameters: 
      <return>       Primary productivity in milligrams Carbon per square meter
                     per day

   !Dependencies:
      function austinPetzold_1986 ( double lambda, double K490 )

         given a reference k490 vlaue, determine k(lambda) for a specified lambda

         ref:
            Austin, R. W., and T. J. Petzold (1986), Spectral dependence of the diffuse
            attenuation coefficient of light in ocean waters, Opt. Eng., 25, 473 – 479

   !Revision History:  

   08-16-2010 first release version (Robert O'Malley)
      [original code written in matlab by T. Westberry]

   01-05-2011   O'Malley
      add uMax trap on mu[m]
      correct z_eu determination

   !References and Credits   
      Westberry, T. Behrenfeld, M.J., Siegel, D.A., and Boss, E.; 2008.  Carbon-based
      primary productivity modeling with vertically resolved photoacclimation.  Global
      Biogeochemical Cycles, Vol. 22, GB2024, doi:10.1029/2007GB003078
      */

  double austinPetzold_1986( double, double );

  double uMax;			/* max growth rate */
  double chlCarbonMax;		/* max chl:carbon ration */
  double nutTempFunc;		/* f(nut,T) */
  double chlCarbonSat;		/* satalite chl:carbon ratio */
  double carbon;		/* bbp converted to carbon */
  double IgFunc;		/* f(Ig) */
  double IgFuncz;               /* f(Ig) below the mixed layer depth */
  double z_eu;			/* euphotic depth at 1% light level */
  double npp;                   /* net primary production */

/* --------------------- */
/*   spectral variables  */
/* --------------------- */

  double lambda[] = { 400, 412, 443, 490, 510, 555, 625, 670, 700 };
  double parFraction[] = { 0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024 };
  double X[] = { .11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000 };
  double e[] = { .64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000 }; 
  double Kw[]= { .01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438 };
  double Kd[9];
  double Kbio;
  double Kdif[9];

  double Klambda[9];
  double Eo[9];
  double Ez_mld[9];
  double par_mld;
  double delChlC;

  double y0;

/* --------------------------- */
/*   depth resolved variables  */
/* --------------------------- */

  double z[200];               /* depths */
  double chl_C[200];           /* chl:c ratio */
  double chlz[200];            /* chl */
  double mu[200];              /* growth */
  double Ezlambda[9][200];     /* fraction of light at nine wavelengths */
  double parz[200];            /* total light */
  double prcnt[200];           /* percent light */
  double Cz[200];              /* carbon */
  double ppz[200];             /* npp */

  int i;
  int m;
  int mzeu;
  double r;
  double prcnt0;
  double prcnt1;
  double z0;
  double z1;
  double numerator;
  double denominator;
  double fraction;
  double deltaZ;

  if(irr <= 0.0){
    return 0.0;
  }

  /* --------------------- */
  /*   initialize values   */
  /* --------------------- */

  z_eu = -9999;     //  1.05.2011
  y0 = 0.0003;                     /* min  chl:c  when  mu = 0 */
  for (i = 0; i<200; i++){
    z[i]= (float)(i+1);
  }
  r = 0.1;
  
  uMax = 2.0;                      /* after Banse (1991) */
  npp = 0.0;
  mzeu = 0;

  for (i=0; i<9; i++) {
    Klambda[i]= austinPetzold_1986(lambda[i],k490);
    Eo[i] = irr * parFraction[i];
    Ez_mld[i] = Eo[i] * 0.975 * exp(-Klambda[i] * mld / 2.0);
  }

  /* ----------------------------- */
  /*   reintegrate to get par at   */
  /*   depth ...                   */
  /*   do trapezoidal integration  */
  /* ----------------------------- */

  par_mld = 0.0;
  for (i = 0; i < 8; i++ ){
    par_mld += (lambda[i+1]-lambda[i])*(Ez_mld[i+1]+Ez_mld[i])/2;
  }

  par_mld /= daylength;

  IgFunc = 1 - exp(-5.0 * par_mld);

  if(bbp < 0.00035)
    bbp = 0.00036;
  carbon = 13000.0 * (bbp - 0.00035);
  
  chlCarbonSat = chl / carbon;

  if ( chlCarbonSat < y0 ) {
    chlCarbonSat = y0;
  }

  chlCarbonMax = 0.022 + (0.045-0.022) * exp(-3.0 * par_mld);
  delChlC = chlCarbonMax - chlCarbonSat;

  nutTempFunc = (chlCarbonSat - y0) / (chlCarbonMax - y0);

  /* ''''''''''''''''''''''''' */
  /*   calculate Kd offset     */
  /*   carry through to depth  */
  /*   non-chl attenuation     */
  /* ------------------------- */

  for (i=0; i<9; i++) {
    Kbio = X[i] * pow(chl,e[i]);
    Kd[i] = Kw[i] + Kbio;
    Kdif[i] = Klambda[i] - Kd[i];
  }

  /* ''''''''''''''''''''''''''''''''''' */
  /*   integrate down the water column   */
  /*   in one-meter steps                */
  /* ----------------------------------- */

  for (m=0; m<200; m++) {

    /* ---------------------------------------------- */
    /*   if you are in the mixed layer, do this way   */
    /* ---------------------------------------------- */

    if ( z[m] < mld ) {
      chl_C[m] = chlCarbonSat;
      chlz[m] = chl_C[m] * carbon;
      mu[m] = uMax * nutTempFunc * IgFunc;

      if ( mu[m] > uMax ) {     //  1.05.2011
	mu[m] = uMax;           //  1.05.2011
      }                         //  1.05.2011
  
      for ( i=0; i<9; i++) {
	Ezlambda[i][m] = Eo[i]*0.975*exp(-Klambda[i]*z[m]);
      }

      parz[m] = 0.0;
      for (i = 0; i < 8; i++ ){
	parz[m] += (lambda[i+1]-lambda[i])*(Ezlambda[i+1][m]+Ezlambda[i][m])/2;
      }

      Cz[m] = carbon;

    } else {

      /* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
      /*   if below mixed layer must treat properties differently   */
      /* ---------------------------------------------------------- */

      for (i=0; i<9; i++) {
	Kbio = X[i] * pow(chlz[m-1],e[i]);     /*  after Morel & Maritorena (2001)  */
	Kd[i] = Kw[i] + Kbio;
	Kd[i] += Kdif[i];
	Ezlambda[i][m] = Ezlambda[i][m-1]*exp(-Kd[i]*1.0);
      }

      parz[m] = 0.0;
      for (i = 0; i < 8; i++ ){
	parz[m] += (lambda[i+1]-lambda[i])*(Ezlambda[i+1][m]+Ezlambda[i][m])/2;
      }

      deltaZ = zno3 - z[m];
      if ( deltaZ < 0 ) {
	deltaZ = 0;
      }

      chl_C[m] = (0.022 + (0.045-0.022) * exp(-3.0 * parz[m] / daylength));
      chl_C[m] -= delChlC * (1-exp(-0.075*deltaZ));

      IgFuncz = 1 - exp(-5.0 * parz[m]/daylength);
      mu[m] = uMax * nutTempFunc * IgFuncz;

      if ( mu[m] > uMax ) {     //  1.05.2011
	mu[m] = uMax;           //  1.05.2011
      }                         //  1.05.2011
  
      if (mu[m-1] >= r ) {
	Cz[m] = carbon;
      } else {
	Cz[m] = carbon * mu[m-1] / r;
      }

      chlz[m] = chl_C[m] * Cz[m];

    }
	   
    prcnt[m] = parz[m] / (irr * 0.975);

    /*  track this to get to the euphotic depth  */

    if ( prcnt[m] >= 0.01 ) {
      mzeu = m;
    } else {

      /* ''''''''''''''''''''''''''' */
      /*   now find 1% light depth   */
      /*   in case the user wants    */
      /*   to use this information   */
      /* --------------------------- */

      if (z_eu == -9999 ) {     // 01.05.11
        prcnt0 = prcnt[mzeu];
        prcnt1 = prcnt[mzeu+1];
        z0 = z[mzeu];
        z1 = z[mzeu+1];
        numerator = prcnt0 - 0.01;
        denominator = prcnt0 - prcnt1;
        fraction = numerator / denominator;
        z_eu = z0 + (z1-z0)*fraction;
      }
    }

    ppz[m] = mu[m] * Cz[m];

  }

  /* ------------------------------- */
  /*   do trapezoidal integration    */
  /*   from m = 0 to m = 200         */
  /* ------------------------------- */

  //  note:  186 m is the euphotic depth for pure water

  if ( mzeu < 186 ) {
    npp = 0;
    for (i = 0; i < 199; i++ ){
      npp += (z[i+1]-z[i])*(ppz[i+1]+ppz[i])/2;
    }
  } else {
    npp = -9999;
  }

  return npp;
}

/* =================================================================  */

double austinPetzold_1986 ( double lambda,
                             double K490 ) {

  double wave[] = { 350, 360, 370, 380, 390, 400,
                    410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
                    510, 520, 530, 540, 550, 560, 570, 580, 590, 600,
                    610, 620, 630, 640, 650, 660, 670, 680, 690, 700 };

  double M[] = { 2.1442, 2.0504, 1.9610, 1.8772, 1.8009, 1.7383,
		 1.7591, 1.6974, 1.6108, 1.5169, 1.4158, 1.3077, 1.1982, 1.0955, 1.0000, 0.9118, 
		 0.8310, 0.7578, 0.6924, 0.6350, 0.5860, 0.5457, 0.5146, 0.4935, 0.4840, 0.4903, 
		 0.5090, 0.5380, 0.6231, 0.7001, 0.7300, 0.7301, 0.7008, 0.6245, 0.4901, 0.2891 };

  double Kdw[] = { 0.0510, 0.0405, 0.0331, 0.0278, 0.0242, 0.0217, 
		   0.0200, 0.0189, 0.0182, 0.0178, 0.0176, 0.0176, 0.0179, 0.0193, 0.0224, 0.0280, 
		   0.0369, 0.0498, 0.0526, 0.0577, 0.0640, 0.0723, 0.0842, 0.1065, 0.1578, 0.2409, 
		   0.2892, 0.3124, 0.3296, 0.3290, 0.3559, 0.4105, 0.4278, 0.4521, 0.5116, 0.6514 };

  double l0;
  double l1;
  double k0;
  double k1;
  double m0;
  double m1;
  double kdiff;
  double mdiff;
  double num;
  double den;
  double frac;
  double Kdw_l;
  double M_l;
  double Kd;

  int ref;
  int i;

  // -- INTERPOLATE TO WAVELENGTH OF INTEREST --  //

  for (i = 1; i < 36; i++) {
    if ( wave[i] >= lambda ) {
      l1 = wave[i];
      k1 = Kdw[i];
      m1 = M[i];
      l0 = wave[i-1];
      k0 = Kdw[i-1];
      m0 = M[i-1];
      break;
    }
  }

  num = lambda - l0;
  den = l1 - l0;
  frac = num / den;

  kdiff = k1 - k0;
  Kdw_l = k0 + frac*kdiff;

  mdiff = m1 - m0;
  M_l = m0 + frac*mdiff;
  

  // -- GET REFERENCE WAVELENGTH (=490 FOR NOW) AND APPLY MODEL -- //

  ref = 14;

  Kd = (M_l/M[ref]) * (K490 - Kdw[ref]) + Kdw_l;

  return Kd;

}
