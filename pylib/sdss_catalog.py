import sys
import os
import string
import numpy
import scipy
import fitsio

## Written by Nao Suzuki

# 2022-07-16 (Sat) DR8 is added
# 2022-07-03 (Sun) Update

class SDSSspec:
      def __init__(self,plate,mjd,fiber):
          self.dr='DR17'
          self.ver='v5_13_2'
          self.plate=int(plate)
          self.mjd=int(mjd)
          self.fiber=int(fiber)
          self.strmjd=str("%05i"%(self.mjd))
          self.strplate=str("%04i"%(self.plate))
          self.strfiber=str("%03i"%(self.fiber))
          self.radeg=0.0
          self.decdeg=0.0
          self.z=0.0
          self.zerr=0.0
          self.zsource='zpipe'
          self.zspzbest=0.0
          self.zspzbesterr=0.0
          self.airmass=1.0
          self.alt=90.0
          self.az=0.0
          self.exptime=0.0
          self.nexp=0
          self.thingid=0

      def read(self):
      #2022-04-22 Fri
          # Define Data Directory
          self.fitsfilename=os.environ['SPECTRO_REDUX']+self.ver+'/'+\
          self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
          # Read FITS header
          h=fitsio.read_header(self.fitsfilename,ext=0)
          self.nspec=h['NAXIS2']
          self.coeff0=h['coeff0']
          self.coeff1=h['coeff1']
          self.npix=h['naxis1']
          self.airmass=h['airmass']
          #if(self.plate <=3509) : self.quality=hdulist[0].header['quality']
          self.exptime=h['exptime']
          self.nexp=h['nexp']/4
          self.idlspec2dver=h['VERSCOMB']
          self.idstart=int(h['coeff0']/0.0001)
          self.idend=self.idstart+self.npix-1

          # Reading FITS Image
          self.fits=fitsio.FITS(self.fitsfilename)
          self.flux=self.fits[0][self.fiber-1,:][0]
          self.ivar=self.fits[1][self.fiber-1,:][0]
          self.mask=self.fits[2][self.fiber-1,:][0]
          self.err=self.ivar
          # Wavelength Definition
          self.wave=10.0**(self.coeff0+self.coeff1*numpy.arange(self.npix))
          self.wave_plus =10.0**(self.coeff0+self.coeff1*(numpy.arange(self.npix)+0.5))
          self.wave_minus=10.0**(self.coeff0+self.coeff1*(numpy.arange(self.npix)-0.5))
          self.dwave=self.wave_plus-self.wave_minus

          self.wid=int(h['coeff0']/0.0001)+numpy.arange(self.npix)

      def write(self):
      # Jun 12, 2012 (Tue) 3pm : Apr 22, 2022 (Fri) 22:52
      # Write out ascii file
          asciidir=os.environ['SPECTRO_REDUX']+self.ver+'/ascii/'+self.strplate
          if(os.path.exists(asciidir)==False):
            os.mkdir(asciidir)

          # Filter Out unnecessary files
          if(self.fiber!=1000):
            asciifilename=asciidir+'/spSpec-'+self.strplate+'-'+self.strmjd+'-'+self.strfiber+'.dat'
          elif(self.fiber==1000):
            asciifilename=asciidir+'/spSpec-'+self.strplate+'-'+self.strmjd+'-000.dat'
          ofile=open(asciifilename,'w')

          #print(asciifilename)
          for i in range(len(self.wid)):
             ofile.write("%10.4f"%(self.wave[i])+\
                          "%15.5e"%(self.flux[i])+\
                          "%15.5e"%(self.ivar[i])+\
                          "%15i"%(self.mask[i])+\
                          "%15i"%(self.wid[i])+str("%2s"%('\n')))
          ofile.close()

      def normalize_at7000(self):
          w6900id=38389
          w7000id=38451
          w7100id=38513

          nflux=self.flux[(w6900id-self.wid[0]):(w7100id-self.wid[0])]
          nivar=self.ivar[(w6900id-self.wid[0]):(w7100id-self.wid[0])]
          nflag=numpy.ones(w7100id-w6900id+1,dtype=numpy.int)
          # Flag Good Data Points
          nflag=numpy.where(nivar>0.0,1,0)
          # Flag S/N > 1 :  numpy.where (condition, yes, no)
          nflag=numpy.where(nflux*numpy.sqrt(nivar) > 1.0,1,0)
          # Extract Good Data Points : extract from nflux with nflag
          nflux_good=numpy.compress(nflag,nflux)

          # Taking Median : 
          if(len(nflux_good) > 5):
             nfactor=numpy.median(nflux_good)
             nfactor_std=numpy.std(nflux_good)
          elif(len(nflux) > 1):
             nfactor=numpy.median(nflux)
             nfactor_std=numpy.std(nflux)
          else:
              nfactor=1.0
              nfactor_std=0.0

          self.flux/=nfactor

class spall():
      def __init__(self):
          self.version='v5_13_2'
          self.dr='DR17'
          self.fitstablename=os.environ['SPECTRO_REDUX']+'spAll-'+self.version+'.fits'
          h=fitsio.read_header(self.fitstablename,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)

      def readstar(self):
          columns=['PLATE','MJD','FIBERID','OBJTYPE','CLASS'=='STAR','SUBCLASS']
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          self.platelist=d['PLATE']
          self.mjdlist=d['MJD']
          self.fiberlist=d['FIBERID']
          self.objtypelist=d['OBJTYPE']
          self.classlist=d['CLASS']
          self.subclasslist=d['SUBCLASS']

      def readtest(self):
          columns=['OBJTYPE','CLASS','SUBCLASS','SPECPRIMARY']
          print('Reading ',self.fitstablename)
          print('Number of Spectra is ',len(self.nspec))
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          self.objtype=d['OBJTYPE']
          self.objclass=d['CLASS']
          self.subclass=d['SUBCLASS']
          self.specprimary=d['SPECPRIMARY']

      def read(self):
          columns=['PLATE','MJD','FIBERID',\
          'RA','DEC',\
          'PLATEQUALITY','ZOFFSET','ZWARNING',\
          'LAMBDA_EFF','XFOCAL','YFOCAL','PLUG_RA','PLUG_DEC',\
          'OBJTYPE','CLASS','SUBCLASS','SPECPRIMARY',\
          'SPEC1_G','SPEC1_R','SPEC1_I',\
          'SPEC2_G','SPEC2_R','SPEC2_I',\
          'PSFMAG','PSFMAGERR','THING_ID',\
          'CMODELMAG','CMODELMAGERR',\
          'FIBER2MAG','FIBER2MAGERR',\
          'SPECTROFLUX','SPECTROFLUX_IVAR',\
          'SPECTROSYNFLUX','SPECTROSYNFLUX_IVAR','EXTINCTION',\
          'AIRMASS','OBJID','Z','Z_ERR',\
          'SN_MEDIAN','SN_MEDIAN_ALL',\
          'ELODIE_FILENAME','ELODIE_OBJECT','ELODIE_SPTYPE',\
          'ELODIE_BV','ELODIE_FEH',\
          'ELODIE_TEFF','ELODIE_LOGG','ELODIE_Z','ELODIE_Z_ERR']
          print('Reading ',self.fitstablename)
          print('Number of Rows is ',len(self.rows))
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          #d=fitsio.read(fitstablename,columns=columns)

          self.platelist=d['PLATE']
          self.mjdlist=d['MJD']
          self.fiberlist=d['FIBERID']

          self.ralist=d['RA']
          self.declist=d['DEC']

          self.platequalitylist=d['PLATEQUALITY']
          self.zoffsetlist=d['ZOFFSET']
          self.zwarninglist=d['ZWARNING']

          self.lambdaefflist=d['LAMBDA_EFF']
          self.xfocallist=d['XFOCAL']
          self.yfocallist=d['YFOCAL']

          self.plug_ralist=d['PLUG_RA']
          self.plug_declist=d['PLUG_DEC']

          self.objtypelist=d['OBJTYPE']
          self.classlist=d['CLASS']
          self.subclasslist=d['SUBCLASS']
          self.specprimarylist=d['SPECPRIMARY']

          self.spec1glist=d['SPEC1_G']
          self.spec1rlist=d['SPEC1_R']
          self.spec1ilist=d['SPEC1_I']

          self.spec2glist=d['SPEC2_G']
          self.spec2rlist=d['SPEC2_R']
          self.spec2ilist=d['SPEC2_I']

          self.psfmaglist   =d['PSFMAG']
          self.psfmagerrlist=d['PSFMAGERR']
          self.thing_idlist =d['THING_ID']

          self.cmodelmaglist   =d['CMODELMAG']
          self.cmodelmagerrlist=d['CMODELMAGERR']

          self.fibermaglist   =d['FIBER2MAG']
          self.fibermagerrlist=d['FIBER2MAGERR']

          # DR8 and DR9 definition
          self.spectrofluxlist=d['SPECTROFLUX']
          self.spectroflux_ivarlist=d['SPECTROFLUX_IVAR']
          # Derive mag from observed spectra
          self.spectromaglist=22.5-2.5*numpy.log10(d['SPECTROFLUX'])
          self.spectromagerrlist=2.5*numpy.log10(1.0+1.0/d['SPECTROFLUX']/numpy.sqrt(d['SPECTROFLUX_IVAR']))

          # Synthetic Flux
          self.spectrosynfluxlist=d['SPECTROSYNFLUX']
          self.spectrosynflux_ivarlist=d['SPECTROSYNFLUX_IVAR']
          # Derive mag from synthesized spectra
          self.spectrosynmaglist=22.5-2.5*numpy.log10(d['SPECTROSYNFLUX'])
          self.spectrosynmagerrlist=2.5*numpy.log10(1.0+1.0/d['SPECTROSYNFLUX']/numpy.sqrt(d['SPECTROSYNFLUX_IVAR']))

          self.extinctionlist=d['EXTINCTION']

          self.airmasslist=d['AIRMASS']
          self.bestobjidlist=d['OBJID']
          self.zspzbestlist=d['Z']
          self.zspzbesterrlist=d['Z_ERR']

          self.sn_medianlist    =d['SN_MEDIAN']
          self.sn_medianalllist =d['SN_MEDIAN_ALL']

          self.elodie_filelist  =d['ELODIE_FILENAME']
          self.elodie_objectlist=d['ELODIE_OBJECT']
          self.elodie_sptypelist=d['ELODIE_SPTYPE']
          self.elodie_bvlist    =d['ELODIE_BV']
          self.elodie_fehlist   =d['ELODIE_FEH']
          self.elodie_tefflist  =d['ELODIE_TEFF']
          self.elodie_logglist  =d['ELODIE_LOGG']
          self.elodie_zlist     =d['ELODIE_Z']
          self.elodie_zerrlist  =d['ELODIE_Z_ERR']
          print('Reading spAll FITS file',self.version,' is done')

# 2022-07-16 (Sat) : rewritten from sdss_data.py
class DR8specobj():
      def __init__(self):
          self.version='dr8'
          self.dr='DR8'
          self.fitstablename=os.environ['SPECTRO_REDUX']+'specObj-dr8.fits'
          h=fitsio.read_header(self.fitstablename,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)

      def read(self):
          columns=['PLATE','MJD','FIBERID',\
          'PLUG_RA','PLUG_DEC',\
          'ZWARNING',\
          'XFOCAL','YFOCAL','PLUG_RA','PLUG_DEC',\
          'CLASS','SUBCLASS','SPECPRIMARY',\
          'SPEC1_G','SPEC1_R','SPEC1_I',\
          'SPEC2_G','SPEC2_R','SPEC2_I',\
          'SPECTROFLUX','SPECTROFLUX_IVAR',\
          'SPECTROSYNFLUX','SPECTROSYNFLUX_IVAR',\
          'SPECOBJID','BESTOBJID',\
          'Z','Z_ERR',\
          'SN_MEDIAN',\
          'ELODIE_FILENAME','ELODIE_OBJECT','ELODIE_SPTYPE',\
          'ELODIE_BV','ELODIE_FEH',\
          'ELODIE_TEFF','ELODIE_LOGG','ELODIE_Z','ELODIE_Z_ERR']

          if(self.dr=='DR8'): 
             columns.append('OBJTYPE')
             columns.append('ORIGOBJID')
             columns.append('TARGETOBJID')
          if(self.dr=='DR9'): 
             columns.append('BOSS_SPECOBJ_ID')
             columns.append('SN_MEDIAN_ALL')

          print('Reading ',self.fitstablename)
          print('Number of Rows is ',len(self.rows))
          # Reading FITS TABLE
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)

          if(self.dr=='DR8'): self.objtypelist=d['OBJTYPE']
          if(self.dr=='DR9'): 
             self.boss_specobjidlist=d['BOSS_SPECOBJ_ID']
             self.sn_medianalllist =d['SN_MEDIAN_ALL']

          self.platelist=d['PLATE']
          self.mjdlist=d['MJD']
          self.fiberlist=d['FIBERID']
          self.zwarninglist=d['ZWARNING']

          self.xfocallist=d['XFOCAL']
          self.yfocallist=d['YFOCAL']

          self.plug_ralist=d['PLUG_RA']
          self.plug_declist=d['PLUG_DEC']

          self.classlist=d['CLASS']
          self.subclasslist=d['SUBCLASS']
          self.specprimarylist=d['SPECPRIMARY']

          self.spec1glist=d['SPEC1_G']
          self.spec1rlist=d['SPEC1_R']
          self.spec1ilist=d['SPEC1_I']

          self.spec2glist=d['SPEC2_G']
          self.spec2rlist=d['SPEC2_R']
          self.spec2ilist=d['SPEC2_I']

          self.spectrofluxlist=d['SPECTROFLUX']
          self.spectroflux_ivarlist=d['SPECTROFLUX_IVAR']
          # Derive mag from observed spectra
          self.specmaglist=22.5-2.5*numpy.log10(d['SPECTROFLUX'])
          self.specmagerrlist=2.5*numpy.log10(1.0+1.0/d['SPECTROFLUX']/numpy.sqrt(d['SPECTROFLUX_IVAR']))

          self.specsynfluxlist=d['SPECTROSYNFLUX']
          self.specsynflux_ivarlist=d['SPECTROSYNFLUX_IVAR']
          # Derive mag from synthesized spectra
          self.specsynmaglist=22.5-2.5*numpy.log10(d['SPECTROSYNFLUX'])
          self.specsynmagerrlist=2.5*numpy.log10(1.0+1.0/d['SPECTROSYNFLUX']/numpy.sqrt(d['SPECTROSYNFLUX_IVAR']))

          self.bestobjidlist=d['BESTOBJID']
          #self.specobjidlist=d['SPECOBJID']
          #self.targetobjidlist=d['TARGETOBJID']
          #self.origobjidlist=d['ORIGOBJID']

          self.zspzbestlist=d['Z']
          self.zspzbesterrlist=d['Z_ERR']
          self.sn_medianlist    =d['SN_MEDIAN']
          self.elodie_filelist  =d['ELODIE_FILENAME']
          self.elodie_objectlist=d['ELODIE_OBJECT']
          self.elodie_sptypelist=d['ELODIE_SPTYPE']
          self.elodie_bvlist    =d['ELODIE_BV']
          self.elodie_fehlist   =d['ELODIE_FEH']
          self.elodie_tefflist  =d['ELODIE_TEFF']
          self.elodie_logglist  =d['ELODIE_LOGG']
          self.elodie_zlist     =d['ELODIE_Z']
          self.elodie_zerrlist  =d['ELODIE_Z_ERR']
          print('Reading specObj FITS file',self.version,' is done')

class DR8photomatchplate():
      def __init__(self):
          self.version='dr8'
          self.dr='DR8'
          self.fitstablename=os.environ['SPECTRO_REDUX']+'photoMatchPlate-dr8.fits'
          h=fitsio.read_header(self.fitstablename,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)

      def read(self):
          columns=['RA','DEC','OBJID','THING_ID']
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          self.ralist=d['RA']
          self.declist=d['DEC']
          self.bestobjidlist=d['OBJID']
          self.thing_idlist=d['THING_ID']
          self.objidlist=d['OBJID']
          del d

# 2022-07-16 (Sat) : rewritten from sdss_data.py
class DR8photoplate():
      def __init__(self):
          self.version='dr8'
          self.dr='DR8'
          self.fitstablename=os.environ['SPECTRO_REDUX']+'photoPlate-dr8.fits'
          h=fitsio.read_header(self.fitstablename,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)

      def read(self):
          columns=['RA','DEC','AIRMASS','PSFMAG','PSFMAGERR','FIBERMAG','FIBERMAGERR',\
                   'CMODELMAG','CMODELMAGERR','THING_ID','OBJID']
          d=fitsio.read(self.fitstablename,columns=columns)
          self.ralist=d['RA']
          self.declist=d['DEC']
          self.airmasslist=d['AIRMASS']
          self.psfmaglist=d['PSFMAG']
          self.psfmagerrlist=d['PSFMAGERR']

          self.fibermaglist=d['FIBERMAG']
          self.fibermagerrlist=d['FIBERMAGERR']

          self.cmodelmaglist   =d['CMODELMAG']
          self.cmodelmagerrlist=d['CMODELMAGERR']
          self.thing_idlist=d['THING_ID']
          self.objidlist=d['OBJID']
          del d
