J/AJ/146/46     Star Formation Rate in nearby galaxies     (Karachentsev+, 2013)
================================================================================
Star formation properties in the Local Volume galaxies via H{alpha} and 
far-ultraviolet fluxes.
    Karachentsev I.D., Kaisina E.I.
   <Astron. J., 146, 46 (2013)>
   =2013AJ....146...46K
================================================================================
ADC_Keywords: Galaxy catalogs ; Galaxies, nearby ; Morphology
Keywords: catalogs - galaxies: dwarf - galaxies: evolution - stars: formation -
          surveys

Abstract:
    A distance-limited sample of 869 objects from the Updated Nearby
    Galaxy Catalog is used to characterize the star formation status of
    the Local Volume population. We present a compiled list of 1217 star
    formation rate (SFR) estimates for 802 galaxies within 11Mpc, derived
    from the H{alpha} imaging surveys and the GALEX far-ultraviolet
    survey. We briefly discuss some basic scaling relations between SFR
    and luminosity, morphology, HI mass, surface brightness, and the
    environment of the galaxies. About 3/4 of our sample consist of dwarf
    galaxies, for which we offer a more refined classification. We note
    that the specific SFR of nearly all luminous and dwarf galaxies does
    not exceed the maximum value: log(SFR/L_K_)=-9.4[yr^-1^]. Most spiral
    and blue dwarf galaxies have enough time to generate their stellar
    mass during the cosmological time, T_0_, with the observed SFRs. They
    dispose of a sufficient amount of gas to support their present SFRs
    over the next T_0_term. We note that only a small fraction of BCD, Im,
    and Ir galaxies (about 1/20) proceed in a mode of vigorous starburst
    activity. In general, the star formation history of spiral and blue
    dwarf galaxies is mainly driven by their internal processes. The
    present SFRs of E, S0, and dSph galaxies typically have 1/30-1/300 of
    their former activity.

Description:
    A distance-limited sample of 869 objects from the Updated Nearby
    Galaxy Catalog (Karachentsev et al., 2013, cat. J/AJ/145/101) is used
    to characterize the star formation status of the Local Volume
    population. We discuss the available observational data on the current
    rate of star formation in the galaxies with distances D<11Mpc, which
    was determined from the integral flux of the galaxy in the emission
    H{alpha} line or from the far-ultraviolet flux, obtained at the GALEX
    orbital telescope.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl   Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80         .   This file
table3.dat     91       802   List of the nearby galaxies with measured SFR
--------------------------------------------------------------------------------

See also:
 J/AJ/145/101   : Updated nearby galaxy catalog (Karachentsev+, 2013)
 J/ApJS/192/6   : A GALEX UV imaging survey of nearby galaxies (Lee+, 2011)
 J/ApJS/173/185 : GALEX ultraviolet atlas of nearby galaxies (Gil de Paz+, 2007)
 J/ApJS/165/307 : Survey for ionization in neutral gas gal.. I. (Meurer+, 2006) 
 J/AJ/128/2170  : H{alpha} imaging of irregular galaxies (Hunter+, 2004)
 J/A+A/414/23   : CCD H{alpha} and R photometry of 334 galaxies (James+, 2004)
 J/ApJS/147/29  : BRHa data of blue compact dwarf galaxies (Gil De Paz+, 2003)

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
  Bytes Format Units     Label   Explanations
--------------------------------------------------------------------------------
   1- 18  A18   ---      Name    Galaxy name
  20- 21  I2    h        RAh     Hour of Right Ascension (J2000)
  23- 24  I2    min      RAm     Minute of Right Ascension (J2000)
  26- 29  F4.1  s        RAs     Second of Right Ascension (J2000)
      31  A1    ---      DE-     Sign of the Declination (J2000)
  32- 33  I2    deg      DEd     Degree of Declination (J2000)
  35- 36  I2    arcmin   DEm     Arcminute of Declination (J2000)
  38- 39  I2    arcsec   DEs     Arcsecond of Declination (J2000)
  41- 42  I2    ---      TT      [-3/11] de Vaucouleurs morphological type (1)
  44- 49  F6.2  mag      BMag    ? Absolute B band magnitude (2)
      51  A1    ---    l_SFRa    [<] Upper limit flag on SFRHa
  52- 56  F5.2 [Msun/yr] SFRa    ? H{alpha} derived integral Star Formation 
                                   Rate (2)
      58  A1    ---    l_Pa      [<] Upper limit flag on PHa
  59- 63  F5.2 [yr-1]    Pa      ? H{alpha} derived evolutionary P parameter (3)
      65  A1    ---    l_Fa      [>] Lower limit flag on FHa
  66- 70  F5.2  [-]      Fa      ? H{alpha} derived evolutionary F parameter (3)
      72  A1    ---    l_SFRu    [<] Upper limit flag on SFRFUV
  73- 77  F5.2 [Msun/yr] SFRu    ? FUV derived integral Star Formation Rate (2)
      79  A1    ---    l_Pu      [<] Upper limit flag on PFUV
  80- 84  F5.2 [yr-1]    Pu      ? FUV derived evolutionary P parameter (3)
      86  A1    ---    l_Fu      [>] Lower limit flag on FFUV
  87- 91  F5.2  [-]      Fu      ? FUV derived evolutionary F parameter (3)
--------------------------------------------------------------------------------
Note (1): From de Vaucouleurs et al. (1991, cat. VII/155). Morphological type 
     are: 2-8=Sa-Sdm; 9=irregular magellanic, Im and Blue Compact Dwarf (BCD); 
     10=dwarf irregular, Ir. See Table1 for more details.
Note (2): Corrected for the Galactic and internal extinction.
Note (3): Corrected for the Galactic and internal extinction. To characterize 
     the evolutionary status of a sample of galaxies, Karachentsev & Kaisin 
     (2007AJ....133.1883K) proposed a diagnostic "past-future" (PF) diagram, 
     where the dimensionless parameters P=log(SFR.T_0_/L_K_) and 
     F=log(1.85M_H1_/SFR.T_0_) are independent from errors in finding
     distances  to the galaxies. The parameter P is actually the sSFR
     (specific star formation rate) over the entire age  scale of the
     universe, T_0_=13.7Gyr. The F parameter corresponds to the  notion
     of gas depletion time, expressed in units of T_0_. The coefficient
     1.85 at M_HI_ is introduced in order to account for the contribution
     of helium and molecular hydrogen to the total mass of the gas
     (Fukugita & Peebles, 2004ApJ...616..643F).
--------------------------------------------------------------------------------

History:
    From electronic version of the journal 

================================================================================
(End)                Greg Schwarz [AAS], Sylvain Guehenneux [CDS]    30-Jun-2014
