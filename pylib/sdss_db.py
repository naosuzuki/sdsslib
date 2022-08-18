import sys
import os
import string
import numpy
import scipy
import fitsio
from astropy.io import fits
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import extinction
import astrolib
import cosmology

## Written by Nao Suzuki

# 2022-08-16 (Tue) deredshift is added
# 2022-08-13 (Sat) Create 2D Fits is added
# 2022-07-16 (Sat) DR8 is added
# 2022-07-03 (Sun) Update

class SDSSspec:
      def __init__(self,plate,mjd,fiber):
          if(int(plate)>3509): 
             self.dr='DR17'
             self.ver='v5_13_2'
          elif(int(plate)<=3509): 
             self.dr='DR8'
             self.ver='DR8'
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
          self.fitsfilename='sdss.fits'

      def deredshift(self,dlog=0.0001):
          # Bring spectra into restframe 
          # Written in 2012, revised 2022-08-17 (Wed)

          # Deredshift
          self.restwave=self.wave/(1.0+self.z)
          # Log Linear Velocity Binning : vpix=70.0 km/s
          self.dlog=dlog

          # Log Liner Rebinning : 
          [self.rwave,self.rflux,self.rivar,self.rid,self.rmask]=\
          astrolib.exec_logrebinning_fluxivarmask(dlog,self.restwave,self.flux,self.ivar,self.mask)

      def read(self):
      #2022-04-22 Fri
          # Define Data Directory
          if(self.dr=='DR17'):
             self.fitsfilename=os.environ['SPECTRO_REDUX']+self.ver+'/'+\
             self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
          elif(self.dr=='DR8'):
             tmpfitsfile=os.environ['SPECTRO_REDUX']+self.ver+'/'+\
             self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
             if(os.path.exists(tmpfitsfile)): 
                self.fitsfilename=tmpfitsfile
             else:
                tmpfitsfile=os.environ['SPECTRO_REDUX']+self.ver+'/103/'+\
                self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
                if(os.path.exists(tmpfitsfile)): 
                   self.fitsfilename=tmpfitsfile
                else: 
                   tmpfitsfile=os.environ['SPECTRO_REDUX']+self.ver+'/104/'+\
                   self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
                   if(os.path.exists(tmpfitsfile)): 
                      self.fitsfilename=tmpfitsfile
                   else:
                      tmpfitsfile=os.environ['SPECTRO_REDUX']+self.ver+'/26/'+\
                      self.strplate+'/spPlate-'+self.strplate+'-'+self.strmjd+'.fits'
                      if(os.path.exists(tmpfitsfile)): 
                         self.fitsfilename=tmpfitsfile
                      else:
                         print('DR8 data does not exist')
                         return
              
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

          # Plate Properties
          self.lambdaeff=5400.
          self.xfocal=0.0
          self.yfocal=0.0
          self.zoffset=0.0

          self.sdssclass='STAR'
          self.subclass=''

      def photopoints20(self):
          # Find 20 photometric points from GAIA XP spectrum
          # We have 251 spectral points = 11pix*1 + 12pix*19 = 20 points
          self.ptx=numpy.zeros(20)
          self.pty=numpy.zeros(20)
          self.ptyerr=numpy.zeros(20)
          self.ptmask=numpy.zeros(20,dtype=numpy.int32)
          self.ptwid=numpy.zeros(20,dtype=numpy.int32)

          # S/N > 5 criterion for good flux value
          spec_mask=numpy.where(self.flux*numpy.sqrt(self.ivar)>5.0,1.0,0)
          # Flux area = flux * dwave
          spec_maskedflux=self.flux*spec_mask*self.dwave
          spec_maskedwave=self.dwave*spec_mask

          wid0=35286        ; wid1=35646
          wid0_end=wid1-120 ; wid1_start=wid1-119 ; wid1_end=wid1+120
          # Calculate the first point
          self.ptx[0]=10.0**(self.coeff1*float(wid0)) ; self.ptwid[0]=wid0

          # Define the first pixel
          if(self.wid[0]<wid0_end):
            nbin=wid0_end-self.wid[0]+1
            if(numpy.sum(spec_mask[0:nbin])>0): 
               self.pty[0]=numpy.sum(spec_maskedflux[0:nbin])/numpy.sum(spec_maskedwave[0:nbin])
               self.ptmask[0]=1

          # Define the pixel 1-19
          for j in range(1,20):
             self.ptwid[j]=wid1+240*(j-1)
             self.ptx[j]=10.0**(self.coeff1*float(self.ptwid[j]))
             wid_start=self.ptwid[j]-119
             wid_end  =self.ptwid[j]+120

             if((self.wid[0] < wid_end) and (self.wid[-1] > wid_start)):
                if(self.wid[0] >= wid_start):
                   nstart=0 ;  nend=wid_end-self.wid[0]+1
                elif((self.wid[0] < wid_start) and (self.wid[-1] >= wid_end)):
                   nstart=wid_start-self.wid[0] ; nend=wid_end-self.wid[0]+1
                elif((self.wid[-1] >= wid_start) and (self.wid[-1] < wid_end)):
                   nstart=wid_start-self.wid[0] ; nend=self.wid[-1]-self.wid[0]+1

                if(numpy.sum(spec_mask[nstart:nend])>0): 
                   self.pty[j]=numpy.sum(spec_maskedflux[nstart:nend])/numpy.sum(spec_maskedwave[nstart:nend])
                   self.ptmask[j]=1
                   self.ptyerr[j]=self.pty[j]/numpy.average(self.flux[nstart:nend]*numpy.sqrt(self.ivar[nstart:nend]))
        
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

def create_2dspec(df,fitsfilename,objtype,flag_gaia,flag_restframe):
# 2022-08-12 LBNL
# Create 2D FITS File from DataFrame

   # Wavelength Range : 3536.71 - 10415.98
   coeff0=3.5486 
   coeff1=0.0001
   startID=35486
   endID  =40177
   npixall=endID-startID+1

   if(flag_restframe==True):
   # Distance Modulus I : Setup
      ndata=len(df)
      mu_list=numpy.zeros(ndata,numpy.float32)
      muerr_list=numpy.zeros(ndata,numpy.float32)
      muz2_list=numpy.zeros(ndata,numpy.float32)
      muz2err_list=numpy.zeros(ndata,numpy.float32)
      startID_list=numpy.zeros(ndata,numpy.int32)
      npix_list=numpy.zeros(ndata,numpy.int32)
      coeff0_list=numpy.zeros(ndata,numpy.float32)

      cosmo=cosmology.Cosmology()
      cosmo.h=0.7
      cosmo.Ol=0.7
      cosmo.Om=0.3
      z1=2.0
      z2mu=cosmo.DistMod(z1)

   #wave=10.0**(coeff0+numpy.arrange(npixall)*coeff1)
   # Define Arrays
   imageflux=numpy.zeros((len(df),npixall),dtype=numpy.float32)
   imageivar=numpy.zeros((len(df),npixall),dtype=numpy.float32)
   imagemask=numpy.zeros((len(df),npixall),dtype=numpy.int32)

   plate_list =df['plate'].to_numpy()
   mjd_list   =df['mjd'].to_numpy()
   fiber_list   =df['fiber'].to_numpy()

   ra_list =df['ra'].to_numpy()
   dec_list=df['dec'].to_numpy()

   thing_id_list=df['thing_id'].to_numpy()
   snall_list   =df['snall'].to_numpy()

   if((objtype=='star') or (objtype=='quasar')):
   # PSF Magnitudes for Point Source
      psfmag_u_list=df['psfmag_u'].to_numpy()
      psfmagerr_u_list=df['psfmagerr_u'].to_numpy()
      psfmag_g_list=df['psfmag_g'].to_numpy()
      psfmagerr_g_list=df['psfmagerr_g'].to_numpy()
      psfmag_r_list=df['psfmag_r'].to_numpy()
      psfmagerr_r_list=df['psfmagerr_r'].to_numpy()
      psfmag_i_list=df['psfmag_i'].to_numpy()
      psfmagerr_i_list=df['psfmagerr_i'].to_numpy()
      psfmag_z_list=df['psfmag_z'].to_numpy()
      psfmagerr_z_list=df['psfmagerr_z'].to_numpy()
   elif(objtype=='galaxy'):
   # Fiber Magnitudes for Galaxies
      fibermag_u_list=df['fibermag_u'].to_numpy()
      fibermagerr_u_list=df['fibermagerr_u'].to_numpy()
      fibermag_g_list=df['fibermag_g'].to_numpy()
      fibermagerr_g_list=df['fibermagerr_g'].to_numpy()
      fibermag_r_list=df['fibermag_r'].to_numpy()
      fibermagerr_r_list=df['fibermagerr_r'].to_numpy()
      fibermag_i_list=df['fibermag_i'].to_numpy()
      fibermagerr_i_list=df['fibermagerr_i'].to_numpy()
      fibermag_z_list=df['fibermag_z'].to_numpy()
      fibermagerr_z_list=df['fibermagerr_z'].to_numpy()

   if(objtype=='star'):
   # Stars
      object_list=df['object'].to_numpy()
      sptype_list=df['sptype'].to_numpy()
      bv_list    =df['bv'].to_numpy()
      feh_list   =df['feh'].to_numpy()
      teff_list  =df['teff'].to_numpy()
      logg_list  =df['logg'].to_numpy()
      if(flag_gaia):
         parallax_list   =df['parallax'].to_numpy()
         parallaxerr_list=df['parallax_error'].to_numpy()
         parallaxsnr_list=df['parallax_over_error'].to_numpy()
         gaia_gmag_list  =df['phot_g_mean_mag'].to_numpy()
         gaia_bpmag_list =df['phot_bp_mean_mag'].to_numpy()
         gaia_rpmag_list =df['phot_rp_mean_mag'].to_numpy()
         gaia_bprp_list  =df['bp_rp'].to_numpy()
         gaia_mg_list    =gaia_gmag_list+5.0*numpy.log10(parallax_list)-10.0
         gaia_pc_list    =1000.0/parallax_list
         # Replace NAN value by -1
         gaia_mg_list    =numpy.nan_to_num(gaia_mg_list,nan=-1.0)
         gaia_pc_list    =numpy.nan_to_num(gaia_pc_list,nan=-1.0)
         #phot_g_mean_mag+5*log10(parallax)-10 as mg, 1000/parallax as dpc
   else:
   # Extragalacic Objects
      z_list       =df['z'].to_numpy()
      zerr_list    =df['zerr'].to_numpy()
      zwarning_list=df['zwarning'].to_numpy()
   # E(B-V) of Milky Way
      coords  =SkyCoord(ra_list,dec_list,frame='icrs',unit='deg')
      sfd     =SFDQuery()
      ebv_list=sfd(coords)

   Rv=3.1
   for i in range(len(df)):
      plate=df['plate'].iloc[i]
      mjd  =df['mjd'].iloc[i]
      fiber=df['fiber'].iloc[i]
   
      # Define Spectrum
      spec=SDSSspec(plate,mjd,fiber)
      spec.read()

      if(objtype!='star'):
         spec.z=z_list[i]
         #Milky Way Extinction Correction by CCM89
         extinction_mag=extinction.ccm89(spec.wave,ebv_list[i]*Rv,Rv,unit='aa')
         newflux=extinction.remove(extinction_mag,spec.flux,inplace=False)
         ivar1=extinction.apply(extinction_mag,spec.ivar,inplace=False)
         ivar2=extinction.apply(extinction_mag,ivar1,inplace=False)
         spec.flux=newflux  ; spec.ivar=ivar2 
         del newflux ; del ivar1 ; del ivar2

         if(flag_restframe==True):
         # Deredshift the Spectrum
            spec.deredshift(dlog=0.0001)
            specmask1=numpy.where((spec.rmask == 0 ),1.0,0.0)
            specmask2=numpy.where((spec.rmask == 2**24 ),1.0,0.0)
            specmask3=specmask1+specmask2
            specmask=numpy.where(specmask3 >=1.0,1.0,0.0)
            spec.rivar*=specmask
            del specmask1 ; del specmask2 ; del specmask3 ; del specmask

            startID_list[i]=spec.rid[0]
            coeff0_list[i]="%6.4f"%(spec.rid[0]*coeff1)
            npix_list[i]=len(spec.rid)

         # Distance Modulus II: Calculations
            z=z_list[i]
            [dmu]=cosmology.find_dmuerr(z)
            [dmuz2]=cosmology.find_dmuerrz2(z)
            mu=cosmo.DistMod(z)
            mu_list[i]=mu
            muerr_list[i]=dmu
            muz2_list[i]=mu-z2mu
            muz2err_list[i]=dmuz2

      # Define the spectrum range : first pixel and the last pixel
      if(spec.wid[0] > startID):
         kstart=0
         jstart=spec.wid[0]-startID
      else:
         kstart=spec.wid[0]-startID
         jstart=0

      if(spec.wid[-1] < endID):
         kend=spec.npix-1
         jend=spec.wid[-1]-startID
      else:
         kend=endID-spec.wid[0]
         jend=npixall-1

      #print('reading',i,jstart,jend)
      print('reading',i,plate,mjd,fiber,'fraction=',"%5.2f"%(i/len(df)))
      if(flag_restframe==False):
         imageflux[i,jstart:jend]=spec.flux[kstart:kend]
         imageivar[i,jstart:jend]=spec.ivar[kstart:kend]
         imagemask[i,jstart:jend]=spec.mask[kstart:kend]
      elif(flag_restframe==True):
         # For Rest Frame Data
         jstart=0 ; jend=len(spec.rid)-1
         kstart=0 ; kend=len(spec.rid)-1
         imageflux[i,jstart:jend]=spec.flux[kstart:kend]
         imageivar[i,jstart:jend]=spec.ivar[kstart:kend]
         imagemask[i,jstart:jend]=spec.mask[kstart:kend]

   hdu1=fits.PrimaryHDU(imageflux)
   hdu2=fits.ImageHDU(imageivar)
   hdu3=fits.ImageHDU(imagemask)

   col1=fits.Column(name='ra',format='E',array=ra_list)
   col2=fits.Column(name='dec',format='E',array=dec_list)
   col3=fits.Column(name='thing_id',format='K',array=thing_id_list)
   col4=fits.Column(name='plate',format='J',array=plate_list)
   col5=fits.Column(name='mjd',format='J',array=mjd_list)
   col6=fits.Column(name='fiber',format='J',array=fiber_list)
   col7=fits.Column(name='snr',format='E',array=snall_list)

   if((objtype=='star') or (objtype=='quasar')):
      col8=fits.Column(name='psfmag_u',format='E',array=psfmag_u_list)
      col9=fits.Column(name='psfmagerr_u',format='E',array=psfmagerr_u_list)
      col10=fits.Column(name='psfmag_g',format='E',array=psfmag_g_list)
      col11=fits.Column(name='psfmagerr_g',format='E',array=psfmagerr_g_list)
      col12=fits.Column(name='psfmag_r',format='E',array=psfmag_r_list)
      col13=fits.Column(name='psfmagerr_r',format='E',array=psfmagerr_r_list)
      col14=fits.Column(name='psfmag_i',format='E',array=psfmag_i_list)
      col15=fits.Column(name='psfmagerr_i',format='E',array=psfmagerr_i_list)
      col16=fits.Column(name='psfmag_z',format='E',array=psfmag_z_list)
      col17=fits.Column(name='psfmagerr_z',format='E',array=psfmagerr_z_list)
   elif(objtype=='galaxy'):
      col8=fits.Column(name='fibermag_u',format='E',array=fibermag_u_list)
      col9=fits.Column(name='fibermagerr_u',format='E',array=fibermagerr_u_list)
      col10=fits.Column(name='fibermag_g',format='E',array=fibermag_g_list)
      col11=fits.Column(name='fibermagerr_g',format='E',array=fibermagerr_g_list)
      col12=fits.Column(name='fibermag_r',format='E',array=fibermag_r_list)
      col13=fits.Column(name='fibermagerr_r',format='E',array=fibermagerr_r_list)
      col14=fits.Column(name='fibermag_i',format='E',array=fibermag_i_list)
      col15=fits.Column(name='fibermagerr_i',format='E',array=fibermagerr_i_list)
      col16=fits.Column(name='fibermag_z',format='E',array=fibermag_z_list)
      col17=fits.Column(name='fibermagerr_z',format='E',array=fibermagerr_z_list)

   if(objtype=='star'):
      col18=fits.Column(name='object',format='A',array=object_list)
      col19=fits.Column(name='sptype',format='A',array=sptype_list)
      col20=fits.Column(name='bv',format='E',array=bv_list)
      col21=fits.Column(name='feh',format='E',array=feh_list)
      col22=fits.Column(name='teff',format='E',array=teff_list)
      col23=fits.Column(name='logg',format='E',array=logg_list)
      cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,\
                        col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,\
                        col21,col22,col23])
      if(flag_gaia):
         col24=fits.Column(name='parallax',format='E',array=parallax_list)
         col25=fits.Column(name='parallaxerr',format='E',array=parallaxerr_list)
         col26=fits.Column(name='parallaxsnr',format='E',array=parallaxsnr_list)
         col27=fits.Column(name='gaia_gmag',format='E',array=gaia_gmag_list)
         col28=fits.Column(name='gaia_bpmag',format='E',array=gaia_bpmag_list)
         col29=fits.Column(name='gaia_rpmag',format='E',array=gaia_rpmag_list)
         col30=fits.Column(name='gaia_bprp',format='E',array=gaia_bprp_list)
         col31=fits.Column(name='gaia_mg',format='E',array=gaia_mg_list)
         col32=fits.Column(name='gaia_pc',format='E',array=gaia_pc_list)
         del cols
         cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,\
              col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,\
              col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,\
              col31,col32])

   else:
      col18=fits.Column(name='z',format='E',array=z_list)
      col19=fits.Column(name='zerr',format='E',array=zerr_list)
      col20=fits.Column(name='zwarning',format='J',array=zwarning_list)
      col21=fits.Column(name='E(B-V)',format='E',array=ebv_list)
      if(flag_restframe==False):
         cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,\
             col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21])
      elif(flag_restframe==True):
         col22=fits.Column(name='COEFF0',format='E',array=coeff0_list)
         col23=fits.Column(name='STARTID',format='J',array=startID_list)
         col24=fits.Column(name='DM',format='E',array=mu_list)
         col25=fits.Column(name='DMERR',format='E',array=muerr_list)
         col26=fits.Column(name='DMZ2',format='E',array=muz2_list)
         col27=fits.Column(name='DMZ2ERR',format='E',array=muz2err_list)
         cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,\
                   col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,\
                   col21,col22,col23,col24,col25,col26,col27])

#  Define FITS Table
   tbhdu=fits.BinTableHDU.from_columns(cols)

#  Define FITS Image and Extention
   hdulist=fits.HDUList([hdu1,hdu2,hdu3,tbhdu])

#  Writing FITS Header
   hdr=hdulist[0].header
   hdr.set('COEFF0',coeff0)
   hdr.set('COEFF1',coeff1)
   hdr['comment']='Created by Nao Suzuki 2022-08-12'
   hdr['comment']='1st EXT=Flux, 2nd EXT=INVVAR, 3rd EXT=MASK'
   hdr['comment']='SDSS DR17'
   hdr['comment']='Wavelength=1.0**(COEFF0+COEFF1*i)'
   hdr['comment']=''
   hdr['comment']='1:  RA          (SDSS)'
   hdr['comment']='2:  DEC         (SDSS)'
   hdr['comment']='3:  THING_ID    (SDSS)'
   hdr['comment']='4:  PLATE       (SDSS)'
   hdr['comment']='5:  MJD         (SDSS)'
   hdr['comment']='6:  FIBER       (SDSS)'
   hdr['comment']='7:  S/N         (SDSS)'
   hdr['comment']='8:  PSFMAG_u    (SDSS)'
   hdr['comment']='9:  PSFMAGerr_u (SDSS)'
   hdr['comment']='10: PSFMAG_g    (SDSS)'
   hdr['comment']='11: PSFMAGerr_g (SDSS)'
   hdr['comment']='12: PSFMAG_r    (SDSS)'
   hdr['comment']='13: PSFMAGerr_r (SDSS)'
   hdr['comment']='14: PSFMAG_i    (SDSS)'
   hdr['comment']='15: PSFMAGerr_i (SDSS)'
   hdr['comment']='16: PSFMAG_z    (SDSS)'
   hdr['comment']='17: PSFMAGerr_z (SDSS)'
   if(objtype=='star'):
      hdr['comment']='18: OBJECT      (SDSS)'
      hdr['comment']='19: SPTYPE      (SDSS)'
      hdr['comment']='20: BV          (SDSS)'
      hdr['comment']='21: FEH         (SDSS)'
      hdr['comment']='22: TEFF        (SDSS)'
      hdr['comment']='23: LOGG        (SDSS)'
      if(flag_gaia):
         hdr['comment']='24: PARALLAX    (GAIA)'
         hdr['comment']='25: PARALLAXERR (GAIA)'
         hdr['comment']='26: PARALLAXSNR (GAIA)'
         hdr['comment']='27: GAIA g-mag  (GAIA)'
         hdr['comment']='28: GAIA bp-mag (GAIA)'
         hdr['comment']='29: GAIA bp-mag (GAIA)'
         hdr['comment']='30: GAIA bp-rp  (GAIA)'
         hdr['comment']='31: GAIA M_g    (GAIA)'
         hdr['comment']='32: GAIA pc     (GAIA)'
   else:
      hdr['comment']='18: Z           (SDSS)'
      hdr['comment']='19: Zerr        (SDSS)'
      hdr['comment']='20: ZWARNING    (SDSS)'
      hdr['comment']='21: E(B-V)      (SDSS)'

   if(os.path.exists(fitsfilename)): os.remove(fitsfilename) 
   hdulist.writeto(fitsfilename)

class SDSSspec2d:
   def __init__(self,fitsfilename):
      # Reading FITS Image
      self.fitsfilename=fitsfilename
      self.hdul=fits.open(fitsfilename)
      hdr =self.hdul[0].header

      self.npix =hdr['NAXIS1']
      self.nspec=hdr['NAXIS2']
      self.coeff0=hdr['COEFF0']
      self.coeff1=hdr['COEFF1']
      self.wID=int(self.coeff0/self.coeff1)+numpy.arange(self.npix)
      self.wave=10.0**(self.wID*self.coeff1)

   def read(self):
      self.imageflux=self.hdul[0].data
      self.imageivar=self.hdul[1].data
      self.imagemask=self.hdul[2].data
      tbl=self.hdul[3].data
      self.ra_list=tbl.field('ra')
      self.dec_list=tbl.field('dec')
      self.thing_id_list=tbl.field('thing_id')
      self.snr_list=tbl.field('snr')
      #self.plate_list=tbl.field('plate')
      #self.mjd_list=tbl.field('mjd')
      #self.fiber_list=tbl.field('fiber')

   def read_gaia(self):
      tbl=self.hdul[3].data
      self.parallax_list=tbl.field('parallax')
      self.parallaxerr_list=tbl.field('parallaxerr')
      self.parallaxsnr_list=tbl.field('parallaxsnr')
      self.gaia_gmag_list=tbl.field('gaia_gmag')
      self.gaia_bpmag_list=tbl.field('gaia_bpmag')
      self.gaia_rpmag_list=tbl.field('gaia_rpmag')
      self.gaia_bprp_list=tbl.field('gaia_bprp')
      self.gaia_mg_list=tbl.field('gaia_mg')
      self.gaia_pc_list=tbl.field('gaia_pc')

   def read_psfmag(self):
      tbl=self.hdul[3].data
      self.psfmag_u_list=tbl.field('psfmag_u')
      self.psfmagerr_u_list=tbl.field('psfmagerr_u')
      self.psfmag_g_list=tbl.field('psfmag_g')
      self.psfmagerr_g_list=tbl.field('psfmagerr_g')
      self.psfmag_r_list=tbl.field('psfmag_r')
      self.psfmagerr_r_list=tbl.field('psfmagerr_r')
      self.psfmag_i_list=tbl.field('psfmag_i')
      self.psfmagerr_i_list=tbl.field('psfmagerr_i')
      self.psfmag_z_list=tbl.field('psfmag_z')
      self.psfmagerr_z_list=tbl.field('psfmagerr_z')

   def read_fibermag(self):
      tbl=self.hdul[3].data
      self.fibermag_u_list=tbl.field('fibermag_u')
      self.fibermagerr_u_list=tbl.field('fibermagerr_u')
      self.fibermag_g_list=tbl.field('fibermag_g')
      self.fibermagerr_g_list=tbl.field('fibermagerr_g')
      self.fibermag_r_list=tbl.field('fibermag_r')
      self.fibermagerr_r_list=tbl.field('fibermagerr_r')
      self.fibermag_i_list=tbl.field('fibermag_i')
      self.fibermagerr_i_list=tbl.field('fibermagerr_i')
      self.fibermag_z_list=tbl.field('fibermag_z')
      self.fibermagerr_z_list=tbl.field('fibermagerr_z')

   def extract(self,i):
   # extract ith spectrum
      self.flux=self.imageflux[i,:]
      self.ivar=self.imageivar[i,:]
      self.mask=self.imagemask[i,:]
