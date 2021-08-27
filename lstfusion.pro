;Initial code for LST

;see https://www.harrisgeospatial.com/docs/image_goes16.html
; see https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf
PRO ch_readgoes
  indir = '/air/incoming/CHEESEHEAD/goes/lst_grids/'
  outdir = '/air/adesai/cheesehead/lstfusion/'
  stub = 'OR_ABI-L2-LSTC-M6_G16_s2019'
                                ;DDDHH
                                ;152 to 304
  doys = numgen(152,304)
  outlat = numgen(46.1,45.8,n=30)
  outlon = numgen(-90.5,-90,n=30)

  lsts_i = make_array(n_elements(doys),24,n_elements(outlon),n_elements(outlat),/float,value=nan())

  FOR d = 0,n_elements(doys)-1 DO BEGIN
    thedoy = doys[d]
    doystr = string(thedoy,format='(i3.3)')
    FOR h = 0,23 DO BEGIN
      hrstr = string(h,format='(i2.2)')
      files = file_search(indir+stub+doystr+hrstr+'*.nc',count=nc)
      IF nc GT 0 THEN BEGIN
        print,'Reading day '+doystr+' Hour '+hrstr
        file = files[0]
        data = NCDF_Parse(file, /READ_DATA)
        lst = data['LST','_DATA']
        dqf = data['DQF','_DATA']
        badlst = where(dqf NE 0,nbl)
        IF nbl GT 0 THEN lst[badlst] = nan()
        lst = screen_arr(lst,200,400)
        IF d EQ 0 AND h EQ 0 THEN BEGIN 
          x = data['x','_DATA']
          y = data['y','_DATA']
          center_lon = data['geospatial_lat_lon_extent','geospatial_lon_nadir', '_DATA']
          proj = map_proj_init('GOES-R',CENTER_LONGITUDE=center_lon)
          lats = make_array(n_elements(x),n_elements(y),/float,value=nan())
          lons = lats
          lsts = make_array(n_elements(doys),24,n_elements(x),n_elements(y),/float,value=nan())
          FOR i = 0,n_elements(x)-1 DO BEGIN
            FOR j = 0,n_elements(y)-1 DO BEGIN
              ll = map_proj_inverse(x[i],y[j],map_structure=proj)
              lats[i,j] = ll[1]
              lons[i,j] = ll[0]
            ENDFOR
          ENDFOR
          lons2=lons
          lats2=lats
          triangulate,lons2,lats2,ltri,sphere=lsph,/degrees
        ENDIF 

        lsts[d,h,*,*] = lst
        IF nbl LT 180 THEN BEGIN
          print,'  Good data ',180-nbl,' interpolating!'
          lst2 = lst
          FOR i = 0,n_elements(lons2[*,0])-1 DO FOR j=0,n_elements(lats2[0,*])-1 DO lst2[i,j] = lst[where(lons EQ lons2[i,j] AND lats EQ lats2[i,j])]
          ltest = griddata(lons2,lats2,lst2,triangles=ltri,/sphere,/grid,xout=outlon,yout=outlat,method='NearestNeighbor',/degrees)
          lsts_i[d,h,*,*] = ltest          
        ENDIF
      ENDIF ELSE BEGIN
        print,'Skipping day '+doystr+' Hour '+hrstr
      ENDELSE 
    ENDFOR 
  ENDFOR 

  save,lsts,lsts_i,lats,lons,outlat,outlon,ltri,lsph,lons2,lats2,doys,filename=outdir+'goes_lst.sav'
  
  
END

PRO ch_getnldas
  doys = expand_arr(numgen(275,304),24)
  doystr = string(doys,format='(i3.3)')
  yyyymmdd = jd_to_dy(doystr,y=2019)

  hhmm = replicate_arr(string(numgen(0,23),format='(i2.2)')+'00',30)
  
  noah = 'https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FNLDAS%2FNLDAS_NOAH0125_H.002%2F2019%2F'+doystr+'%2FNLDAS_NOAH0125_H.A'+yyyymmdd+'.'+hhmm+'.002.grb&FORMAT=bmM0Lw&BBOX=25%2C-125%2C53%2C-67&LABEL=NLDAS_NOAH0125_H.A'+yyyymmdd+'.'+hhmm+'.002.grb.SUB.nc4&SHORTNAME=NLDAS_NOAH0125_H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=002'

  mosaic = 'https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FNLDAS%2FNLDAS_MOS0125_H.002%2F2019%2F'+doystr+'%2FNLDAS_MOS0125_H.A'+yyyymmdd+'.'+hhmm+'.002.grb&FORMAT=bmM0Lw&BBOX=25%2C-125%2C53%2C-67&LABEL=NLDAS_MOS0125_H.A'+yyyymmdd+'.'+hhmm+'.002.grb.SUB.nc4&SHORTNAME=NLDAS_MOS0125_H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=002'

  vic = 'https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FNLDAS%2FNLDAS_VIC0125_H.002%2F2019%2F'+doystr+'%2FNLDAS_VIC0125_H.A'+yyyymmdd+'.'+hhmm+'.002.grb&FORMAT=bmM0Lw&BBOX=25%2C-125%2C53%2C-67&LABEL=NLDAS_VIC0125_H.A'+yyyymmdd+'.'+hhmm+'.002.grb.SUB.nc4&SHORTNAME=NLDAS_VIC0125_H&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=002'

  openw,fl,'/air/incoming/CHEESEHEAD/nldas/url_oct.txt',/get_lun
  printf,fl,noah,format='(a0)'
  printf,fl,mosaic,format='(a0)'
  printf,fl,vic,format='(a0)'
  free_lun,fl
  
END


;step 2
;read nldas, interpolate to 30x30, single array

PRO ch_nldas,doys=doys,hrs=hrs,outlat=outlat,outlon=outlon,output=output,nosave=nosave
  outdir = '/air/adesai/cheesehead/lstfusion/'
  indir = '/air/incoming/CHEESEHEAD/nldas/'
  models = ['mosaic/','noah/','vic/']
  modstr = ['MOS','NOAH','VIC']

                                ;NLDAS_VIC0125_H.A20190711.0800.002.grb.SUB.nc4

  IF n_elements(doys) EQ 0 THEN doys = numgen(152,304)
  IF n_elements(hrs) EQ 0 THEN hrs = findgen(24)
  IF n_elements(outlat) EQ 0 THEN outlat = numgen(46.1,45.8,n=30)
  IF n_elements(outlon) EQ 0 THEN outlon = numgen(-90.5,-90,n=30)

  avsft_i = make_array(n_elements(doys),n_elementS(hrs),n_elements(outlon),n_elements(outlat),3,/float,value=nan())

  
  FOR d = 0,n_elements(doys)-1 DO BEGIN
    thedoy = doys[d]
    doystr = jd_to_dy(thedoy,y=2019)
    dstr = string(thedoy,format='(i3.3)')
    print,'NLDAS day '+doystr
    FOR hh = 0,n_elements(hrs)-1 DO BEGIN
      h = hrs[hh]
      hrstr = string(h,format='(i2.2)')
      FOR m = 0,2 DO BEGIN
        IF thedoy GE 275 AND m EQ 2 THEN BEGIN
          fname = indir+models[m]+'NLDAS_'+modstr[m]+'0125_H.A'+dstr+'.'+hrstr+'00.002.grb.SUB.nc4'          
        ENDIF ELSE BEGIN 
          fname = indir+models[m]+'NLDAS_'+modstr[m]+'0125_H.A'+doystr+'.'+hrstr+'00.002.grb.SUB.nc4'
        ENDELSE 
        files = file_search(fname,count=nc)
        IF nc GT 0 THEN BEGIN
          file = files[0]
          nc = ncdf_open(file)
          IF d EQ 0 AND hh EQ 0 THEN BEGIN 
            ncdf_varget,nc,'lat',lat
            ncdf_varget,nc,'lon',lon
          ENDIF 
          ncdf_varget,nc,'AVSFT',avsft
          ncdf_close,nc
          avsft = screen_arr(avsft,200,400)
          FOR i = 0,n_elements(outlon)-1 DO FOR j = 0,n_elements(outlat)-1 DO avsft_i[d,hh,i,j,m] = avsft[closest(lon,outlon[i]),closest(lat,outlat[j])]
                                ;lon, lat, avsft, radt
                                ;others: nswrs, nlwrs, lhtfl, shtfl,
                                ;gflux, dswrf, dlwrf, albdo, tsoil,
                                ;soilm, rzsm, trans, lai,           
 
        ENDIF ELSE BEGIN
          print,'NO file found'
          stop
        ENDELSE 
      ENDFOR
    ENDFOR
  ENDFOR

  output = reform(avsft_i)

  IF ~keyword_set(nosave) THEN save,avsft_i,filename=outdir+'nldas_lst.sav'
;stop  
  ;calc mean and stddev
;linfit with measure_error? 
                                ;mosaic noah vic

                                ;mosaic: NLDAS_MOS0125_H.A20190601.0000.002.grb.SUB.nc4
                                ;time lon lat AVSFT = LST
                                ;starts 6/1 ends 10/1
                                ;doys 152 to 304

                                ;need to get Oct someday

;vic NLDAS_VIC0125_H.A20190711.1600.002.grb.SUB.nc4
;noah NLDAS_NOAH0125_H.A20190711.1600.002.grb.SUB.nc4

;for each day, read all hours, get AVSFT for CHEAD grids
                                ;regrid 30x30
                                ;save


  
END

PRO ch_goesgapfill
;read goes and nldas grids
  ;for each pixel, do a linear fit to average + stddev (as error)

  outdir = '/air/adesai/cheesehead/lstfusion/'
  restore,outdir+'goes_lst.sav'
  restore,outdir+'nldas_lst.sav'
  
  avsft_m = total(avsft_i,5,/nan)/total(finite(avsft_i),5)
  avsft_s = stddev(avsft_i,dim=5)

  lsts_g = lsts_i

  nlfix = lsts_i
  nlfix[*] = nan()
  nlb = reform(nlfix[0,*,*,*])
  nlm = nlb

  FOR i = 0,29 DO BEGIN
    print,'Column ',i
    FOR j = 0,29 DO BEGIN
      FOR k = 0,23 DO BEGIN 
        goes = lsts_i[*,k,i,j]
        nldas = avsft_m[*,k,i,j]
        nldas_s = avsft_s[*,k,i,j]
        gv = where(finite(goes))
        bias = mean(goes[gv]-nldas[gv])
        nldas += bias
        fitexy,nldas[gv],goes[gv],b,m,x_sigma=nldas_s[gv],y_sigma=1.0
        nldas_fix = m * nldas + b
        lsts_g[*,k,i,j] = merge_array(goes,nldas_fix)
        nlfix[*,k,i,j] = nldas_fix
        nlm[k,i,j] = m
        nlb[k,i,j] = b
      ENDFOR
    ENDFOR
  ENDFOR 

  save,nlfix,nlm,nlb,filename=outdir+'goesfuse_lst_stats.sav'
  stop
  save,lsts_g,filename=outdir+'goesfuse_lst.sav'
stop
  k = reform(transpose(lsts_g,[1,0,2,3]),3672,30,30)

  set_plot,'Z'
  !order=1
  !x.style=1
  !y.style=1
  !p.background=255
  !p.color=0

  doystr = string(doys,format='(i3.3)')
  yyyymmdd = expand_arr(jd_to_dy(doystr,y=2019),24)
  hhmm = replicate_arr(string(numgen(0,23),format='(i2.2)')+'00',153)
  
  for i = 0,3671 do begin
    tvplot,k[i,*,*],zrange=[250,310],xrange=[outlon[0],outlon[29]],yrange=[outlat[29],outlat[0]],title='Date: '+yyyymmdd[i]+' Hour: '+hhmm[i]
;add locations of towers (xyouts)
    o = tvrd()
    l[i,*,*]=o
  endfor
  write_mp4,reverse(l[*,*,*],3),'/home/adesai/test5.mp4',fps=24
  set_plot,'X'
  
  stop
  
END


PRO ch_ecostressproc

                                ;for each ecostress image, read file and cloud mask
                                ;process cloud mask
                                ;nearest neighbor resample to UTM
                                ;apply mask
  ;write out 1000x1000

  
;ECOSTRESS tie point
;-90.500337233322497, 46.100231255142141
;0.00062977741696699997, 0.00062977741696699997
;794 by 477

  edir = '/air/incoming/CHEESEHEAD/ecostress/lst/'
  cdir = '/air/incoming/CHEESEHEAD/ecostress/cloud/'
  infile = 'ECO2LSTE.001_SDS_LST_doy2019267004844_aid0001.tif'
  qcfile = 'ECO2LSTE.001_SDS_QC_doy2019267004844_aid0001.tif'
  infile2 = 'ECO2LSTE.001_SDS_LST_doy2019267053924_aid0001.tif'
;  clfile = 'ECO2CLD.001_SDS_CloudMask_doy2019267004844_aid0001.tif'


  openw,lf,'/air/adesai/cheesehead/lstfusion/inputs/ecostreslog.txt',/get_lun

  clf = file_search(cdir+'ECO2CLD.001_SDS_CloudMask_doy*.tif',count=nc)
  FOR i = 0,nc-1 DO BEGIN
    clfile = clf[i]
    doystr = (reverse(strsplit(clf[i],'_',/extract)))[1]
    efile = edir+'ECO2LSTE.001_SDS_LST_'+doystr+'_aid0001.tif'
    IF file_test(efile,/read) THEN BEGIN
      print,'Found match for ',doystr
      eco = read_tiff(efile,geotiff=gtf)
      eco*=0.02
      eco = screen_Arr(eco,200,400)
      geco = where(finite(eco),ngeco)
      teco = ngeco/float(n_elements(eco))
      IF teco GT 0.66 THEN BEGIN 
      
        eco_cl = read_tiff(clfile,geotiff=gtf2)
        clbit = bytarr(n_elements(eco_cl[*,0]),n_elements(eco_Cl[0,*]))
       ; clbit1 = clbit
        clbit7 = clbit
     
        FOR cl = 0,n_elements(eco_cl[*])-1 DO BEGIN
          IF eco_cl[cl] NE -999 THEN BEGIN 
            clbits = bytebin(eco_cl[cl])
            clbit[cl] = clbits[5] ;final cl (thisone)
          ENDIF ELSE clbit7[clbit] = 1 ;missing
        ENDFOR

        percloud = total(clbit)/float(n_elements(clbit)-total(clbit7))
        
        print,'  Percent cloud is: ',percloud*100.
        
        IF percloud LT 0.5 THEN BEGIN
          outlat = numgen(46.1,45.8,n=600)
          outlon = numgen(-90.5,-90,n=600)
          ecoout = make_Array(n_elements(outlon),n_elements(outlat),/float,value=nan())
          ecoout_cl = make_Array(n_elements(outlon),n_elements(outlat),/byte,value=0)
          
          eco_GT = gtf.modeltiepointtag[3:4]
          eco_cl_GT = gtf2.modeltiepointtag[3:4]
          eco_sl = gtf.modelpixelscaletag[0:1]
          eco_cl_sl = gtf2.modelpixelscaletag[0:1]
          eco_X = n_elements(eco[*,0])
          eco_y = n_elements(eco[0,*])
          eco_cl_X = n_elements(clbit[*,0])
          eco_cl_y = n_elements(clbit[0,*])
          
          inlon = eco_GT[0]+eco_sl[0]*findgen(eco_x)
          inlat = eco_GT[1]-eco_sl[1]*findgen(eco_y)
          
          FOR y = 0,n_elements(outlat)-1 DO BEGIN
            FOR x = 0,n_elements(outlon)-1 DO BEGIN
              ecoout[x,y] = eco[closest(inlon,outlon[x]),closest(inlat,outlat[y])]
              ecoout_cl[x,y] = clbit[closest(inlon,outlon[x]),closest(inlat,outlat[y])]
            ENDFOR
          ENDFOR
          ecoout[where(ecoout_cl)]=nan()
          ecoout *= 10
          beco = where(~finite(ecoout),nbe)
          totgood = string((1.-nbe/360000.)*100,format='(f0.1)')
          IF nbe GT 0 THEN ecoout[beco]=-9999
          ecoout = fix(ecoout)
          outfile = '/air/adesai/cheesehead/lstfusion/inputs/ecostress_'+strmid(doystr,3,4)+'_'+strmid(doystr,7,3)+'_'+strmid(doystr,10,2)+'.dat'
          print,'Writing ',outfile
          openw,fl,outfile,/get_lun
          writeu,fl,ecoout
          free_lun,fl
          printf,lf,'ecostress_'+strmid(doystr,3,4)+'_'+strmid(doystr,7,3)+'_'+strmid(doystr,10,2)+'.dat goodval: '+totgood+'%'
        ENDIF
      ENDIF 
      
    ENDIF ELSE BEGIN
      print,'No match found, sad'
    ENDELSE 
  ENDFOR
  free_lun,lf
  

END

  
FUNCTION ch_readfusion,lat,lon,sday=sday,eday=eday,goes=goes,doys=doys,hrs=hrs,transect=transect,version=version

;for each lat lon, find closest, read time series from sday to eday

;if arg_present(goes) also read in that one

                                ;point_lun for each lat/lon

  IF n_elements(sday) EQ 0 THEN sday = 152
  IF n_elements(eday) EQ 0 THEN eday = 304

  IF n_elements(version) EQ 0 THEN vdir = '' ELSE vdir = version + '/'
  
  ndays = eday-sday+1
  ntms = ndays * 24

  IF ~keyword_set(transect) THEN BEGIN 
    doys = expand_arr(numgen(sday,eday),24)
    hrs = replicate_arr(findgen(24),ndays)
  ENDIF ELSE BEGIN
    npts = min([n_elements(lat),n_elements(lon),n_elements(doys),n_elements(hrs)])
  ENDELSE 
  
  IF n_elements(lat) NE n_elements(lon) THEN BEGIN
    print,'Lat and lon not same number of points'
  ENDIF 


  outlat = numgen(46.1,45.8,n=600)
  outlon = numgen(-90.5,-90,n=600)

  clat = closest(outlat,lat)
  clon = closest(outlon,lon)

  locs = (clon*600+clat)*2
  npts = n_elements(locs)

  IF arg_present(goes) THEN BEGIN
;read goes and extract here
    glat = numgen(46.1,45.8,n=30)
    glon = numgen(-90.5,-90,n=30)
    glat1 = closest(glat,lat)
    glon1 = closest(glon,lon)
    gdir = '/air/adesai/cheesehead/lstfusion/'
    restore,gdir+'goesfuse_lst.sav'
    lsts = reform(transpose(lsts_g,[1,0,2,3]),3672,30,30)
    IF keyword_set(transect) THEN BEGIN
      glocs = (doys-152)*24l + fix(hrs)
      lsts = lsts[glocs,*,*]
      goes = make_Array(npts,/float,value=nan())
      FOR l = 0,npts-1 DO goes[l] = lsts[l,glon1[l],glat1[l]]
    ENDIF ELSE BEGIN 
      gloc1 = (sday-152)*24l
      gloc2 = (eday-152)*24l+23
      lsts = lsts[gloc1:gloc2,*,*]
      goes = make_array(npts,ntms,/float,value=nan())
      FOR l = 0,npts-1 DO goes[l,*] = lsts[*,glon1[l],glat1[l]]
    ENDELSE 
  ENDIF


  IF keyword_set(transect) THEN outarr = make_array(npts,/float,value=nan()) ELSE outarr = make_Array(npts,ntms,/float,value=nan())

  fdir = '/air/adesai/cheesehead/lstfusion/outputs/'+vdir+'fuse_2019_' ;DDD_HH.dat

  fname = ''
  get_lun,fl
  IF keyword_set(transect) THEN BEGIN
    FOR l = 0,npts-1 DO BEGIN
      doy = string(doys[l],format='(i3.3)')
      hh = string(hrs[l],format='(i2.2)')
      oldfname = fname
      fname = fdir + doy + '_' + hh + '.dat'
      IF fname NE oldfname THEN BEGIN
        free_lun,fl
        openr,fl,fname,/get_lun
        oldfname = fname
      ENDIF
      point_lun,fl,locs[l]
      dum = fix(0)
      readu,fl,dum
      outarr[l] = dum
    ENDFOR
    free_lun,fl
  ENDIF ELSE BEGIN   
    FOR d = sday,eday DO BEGIN
      doy = string(d,format='(i3.3)')
      FOR h = 0,23 DO BEGIN
        outloc = (d-sday)*24l+h
        hh = string(h,format='(i2.2)')
        fname = fdir + doy + '_' + hh + '.dat'
        openr,fl,fname,/get_lun
        FOR l = 0,npts-1 DO BEGIN
          point_lun,fl,locs[l]
          dum = fix(0)
          readu,fl,dum
          outarr[l,outloc] = dum
        ENDFOR
        free_lun,fl
      ENDFOR 
    ENDFOR
  ENDELSE 

  bv = where(outarr EQ -9999,nbv)
  IF nbv GT 0 THEN outarr[bv] = nan()
  outarr/=10.
  
  return,outarr
  
END


PRO ch_reproject_coord

  glat = numgen(46.1,45.8,n=600)
  glon = numgen(-90.5,-90,n=600)

  geast = numgen(694300,731600,50)
  gnorth = reverse(numgen(5075900,5109300,50))

  glons = replicate_arr(glon,n_elements(glat))
  glats = expand_arr(glat,n_elements(glon))

  geasts = replicate_arr(geast,n_elements(gnorth))
  gnorths = expand_arr(gnorth,n_elements(geast))

  locs = make_array(n_elements(geast),n_elements(gnorth),/long,value=-1)

  PI = 3.14159265D
  FOURTHPI = PI / 4
  deg2rad = PI / 180.0D
  rad2deg = 180.0D / PI
  ReferenceEllipsoid = 'wgs84'

  k0 = 0.9996D
  a = (ellip(ReferenceEllipsoid))[0] ;EquatorialRadius
  eccSquared = (ellip(ReferenceEllipsoid))[1] ;eccentricitySquared
  e1 = ((1.0D)-sqrt((1.0D)-eccSquared))/((1.0D)+sqrt((1.0D)-eccSquared))

  x = double(geasts) - (500000.0D) ;remove 500,000 meter offset for longitude
  y = double(gnorths)
  UTMZone = '15N'
  UZ = strupcase(strtrim(strcompress(UTMZone,/remove_all),2))
  ZoneNumber = long(UZ)
  ZoneLetter = strmid(UZ,strlen(UZ)-1,1)
  NorthernHemisphere = 1

  LongOrigin = (double(ZoneNumber) - (1.0D))*(6.0D) - (180.0D) + (3.0D) ;+3 puts origin in middle of zone
  eccPrimeSquared = (eccSquared)/((1.0D)-eccSquared)

  M = y / k0
  mu = M/(a*((1.0D)-eccSquared/(4.0D)-(3.0D)*eccSquared*eccSquared/(64.0D)-(5.0D)*eccSquared*eccSquared*eccSquared/(256.0D)))
  phi1Rad = mu+((3.0D)*e1/(2.0D)-(27.0D)*e1*e1*e1/(32.0D))*sin((2.0D)*mu)+((21.0D)*e1*e1/(16.0D)-(55.0D)*e1*e1*e1*e1/(32.0D))*sin((4.0D)*mu)+((151.0D)*e1*e1*e1/(96.0D))*sin((6.0D)*mu)
  phi1 = phi1Rad*rad2deg
  
  N1 = a/sqrt((1.0D)-eccSquared*sin(phi1Rad)*sin(phi1Rad))
  T1 = tan(phi1Rad)*tan(phi1Rad)
  C1 = eccPrimeSquared*cos(phi1Rad)*cos(phi1Rad)
  R1 = a*((1.0D)-eccSquared)/  (((1.0D)-eccSquared*sin(phi1Rad)*sin(phi1Rad))^(1.5D))
  D = x/(N1*k0)

  Lat = phi1Rad-(N1*tan(phi1Rad)/R1)*(D*D/(2.0D)-((5.0D)+(3.0D)*T1+(10.0D)*C1-(4.0D)*C1*C1-(9.0D)*eccPrimeSquared)*D*D*D*D/(24.0D)+((61.0D)+(90.0D)*T1+(298.0D)*C1+(45.0D)*T1*T1-(252.0D)*eccPrimeSquared-(3.0D)*C1*C1)*D*D*D*D*D*D/(720.0D))
  Lat = Lat * rad2deg
  
  Long = (D-((1.0D)+(2.0D)*T1+C1)*D*D*D/(6.0D)+((5.0D)-(2.0D)*C1+(28.0D)*T1-(3.0D)*C1*C1+(8.0D)*eccPrimeSquared+(24.0D)*T1*T1)*D*D*D*D*D/(120.0D))/cos(phi1Rad)
  Long = LongOrigin + Long * rad2deg
  
  ;create geast/gnorth grid
  FOR l = 0,n_elements(locs)-1 DO BEGIN
    clox = closest(glon,long[l])
    cloy = closest(glat,lat[l])
    locs[l] = clox+cloy*n_elements(glon)
  ENDFOR 

  save,geast,gnorth,geasts,gnorths,glat,glon,glats,glons,locs,filename='/air/adesai/cheesehead/lstfusion/coords.sav'
  
END

PRO ch_reproject_write,doy,hr,ascii=ascii,version=version
;read in day/hour
;write out reprojects
  restore,'/air/adesai/cheesehead/lstfusion/coords.sav'

  IF n_elements(version) EQ 0 THEN vdir = '' ELSE vdir = version + '/'
   
  fdir = '/air/adesai/cheesehead/lstfusion/outputs/'+vdir+'fuse_2019_'+string(doy,format='(i3.3)')+'_'+string(hr,format='(i2.2)')+'.dat'
  odir = '/air/adesai/cheesehead/lstfusion/outputs/'+vdir+'utm/fuse_2019_'+string(doy,format='(i3.3)')+'_'+string(hr,format='(i2.2)')
  IF keyword_set(ascii) THEN odir+='.asc' ELSE odir+='.tif'

  input = intarr(n_elements(glon),n_elements(glat))
  openr,fl,fdir,/get_lun
  readu,fl,input
  freE_lun,fl
  input = float(input)
  bv = where(input EQ -9999,nbv)
  IF nbv GT 0 THEN input[bv]=nan()
  input/=10.

  output = fltarr(n_elements(geast),n_elements(gnorth))
  output[*] = input[locs]

  IF keyword_set(ascii) THEN BEGIN
    output*=10
    bv = where(~finite(output),nbv)
    IF nbv GT 0 THEN output[bv]=-9999
    output=fix(output)
    openw,fl,odir,/get_lun
    printf,fl,'NCOLS',n_elements(geast),format='(a-15,i0)'
    printf,fl,'NROWS',n_elements(gnorth),format='(a-15,i0)'
    printf,fl,'XLLCORNER',geast[0],format='(a-15,i0)'
    printf,fl,'YLLCORNER',gnorth[n_elements(gnorth)-1],format='(a-15,i0)'
    printf,fl,'CELLSIZE',50,format='(a-15,i0)'
    printf,fl,'NODATA_VALUE',-9999,format='(a-15,i0)'
    FOR i = 0,n_elements(gnorth)-1 DO printf,fl,output[*,i],format='('+string(n_elements(geast),format='(i0)')+'(i0," "))'
    free_lun,fl
  ENDIF ELSE BEGIN
    geo = {ModelTiepointTag : [0,0, 0, geast[0],gnorth[0], 0],  ModelPixelScaleTag  : [50,50,0], $
           GTModelTypeGeoKey : 1, GTRasterTypeGeoKey : 1, ProjectedCSTypeGeoKey : 32615, $
           GEOGRAPHICTYPEGEOKEY : 4326, PCSCitationGeoKey : 'UTM Zone 15 N with WGS84'}
                               
    write_Tiff,odir,output,/float,geotiff=geo,compress=1

  ENDELSE 
  

;ncols 4

  

END


PRO ch_dronecomp,version=version,lidar=lidar,hyspex=hyspex,jfix=jfix
;read in drone, read in UTM grid
  drone_dr = '/bog/incoming/cheesehead/drone-lst/'
  flights = ['07112019','07112019','07112019',$
             '07122019','07122019','07122019','07122019',$
             '08192019','08192019',$
             '08212019','08212019','08212019','08212019','08212019']
  fnumber = ['1','2','4','1','2','3','4','1','2','1','2','3','4','5']
  doys = [192,192,192,193,193,193,193,231,231,233,233,233,233,233]
  tms = [21,21.5,22,15,15.5,16,16.5,17,21,14,16,17,19,20.5]

  IF n_elements(version) EQ 0 THEN version = 'v20200803'
  IF version EQ '' THEN vdir = '' ELSE vdir = version + '/'

  IF keyword_set(lidar) THEN BEGIN
;    restore,'/air/adesai/cheesehead/lidar/can.sav'
;    li_ux = 696329.300
;    li_uy = 5106591.070

    li_dr='/air/incoming/CHEESEHEAD/PriceCounty_LiDAR/processed_tzheng_20191227/'
    dem_f = 'CH_DEM_30km_1m_clipped'
    dsm_f = 'CH_DSM_30km_1m_clipped'
    dem_c1 = [696329.300, 5106591.070]
    dsm_c1 = [696328.474, 5106591.070]
    openu,dem_fl,li_dr+dem_f,/get_lun
    openu,dsm_fl,li_dr+dsm_f,/get_lun
  ENDIF 
    
  FOR i = 6,n_elements(flights)-1 DO BEGIN
    dr = read_tiff(drone_dr+flights[i]+'_flight'+fnumber[i]+'.tif',geotiff=gdr)
    bv = where(dr LT -200,nbv)
    IF nbv GT 0 THEN dr[bv]=nan()
    dr+=273.15
    doy = doys[i]
    hr = tms[i]
    hrf = fix(hr)
    fname = '/air/adesai/cheesehead/lstfusion/outputs/'+vdir+'utm/fuse_2019_'+string(doy,format='(i3.3)')+'_'+string(hr,format='(i2.2)')+'.tif'
    IF ~file_test(fname,/read) THEN ch_reproject_write,doy,hr,version=version
    IF file_test(fname,/read) THEN fuse = read_tiff(fname,geotiff=gfu) ELSE stop
    IF hr NE hrf THEN BEGIN
      hr2 = hrf+1
      fname2 = '/air/adesai/cheesehead/lstfusion/outputs/'+vdir+'utm/fuse_2019_'+string(doy,format='(i3.3)')+'_'+string(hr2,format='(i2.2)')+'.tif'
      IF ~file_test(fname2,/read) THEN ch_reproject_write,doy,hr2,version=version
      IF file_test(fname2,/read) THEN fuse2 = read_tiff(fname2,geotiff=gfu2) ELSE stop
      fuse1 = fuse
      fuse = (fuse1+fuse2)/2.
    ENDIF 

;extract fuse to drone projection

    dr_sc = gdr.modelpixelscaletag[0]
    fu_sc = gfu.modelpixelscaletag[0]
    dr_c1 = gdr.modeltiepointtag[3:4]
    fu_c1 = gfu.modeltiepointtag[3:4]
    fu_loc_x = fix((dr_c1[0]-fu_c1[0])/fu_sc)
    fu_loc_y = fix((fu_c1[1]-dr_c1[1])/fu_sc)
    fu_span_x = round(fu_loc_x + (n_elements(dr[*,0])*dr_sc/fu_sc))
    fu_span_y = round(fu_loc_y + (n_elements(dr[0,*])*dr_sc/fu_sc))
    fu_ext = fuse[fu_loc_x:fu_span_x,fu_loc_y:fu_span_y]

    dr_fu = congrid(fu_ext,n_elements(dr[*,0]),n_elements(dr[0,*]),/center)
    dr_fu[where(~finite(dr))]=nan()

    IF n_elements(hyspex) NE 0 THEN BEGIN
;if july 13, cut off top 49 pixels
      IF keyword_set(jfix) THEN BEGIN
        dr = dr[*,49:*]
        dr_c1[1]-=dr_sc*49
      ENDIF 
      
      ch_hyspex_find,hyspex,[dr_c1[0],dr_c1[0]+dr_sc*n_elements(dr[*,0])],[dr_c1[1],dr_c1[1]-dr_sc*n_elements(dr[0,*])],locs=hlocs,files=hfl
      stop
      IF hlocs[0,0] GT 0 AND hlocs[0,1] GT 0 THEN BEGIN
        hy = ch_hyspex_rd(hyspex,box=hlocs[*],bands=indgen(474))
        bhy = where(hy LT 0,nbhy)
        IF nbhy GT 0 THEN hy[bhy] = nan()

        ;10 m avg
        hy_dr10 = congrid(hy,n_elements(dr[*,0])/23,n_elements(dr[0,*])/23,474,/center,/interp)
        dr10 = congrid(dr,n_elements(dr[*,0])/23,n_elements(dr[0,*])/23,/interp,/center)
        hy_dr10[where(hy_dr10 EQ 0)]=nan()
        fu10 = congrid(dr_fu,n_elements(dr[*,0])/23,n_elements(dr[0,*])/23,/center)
        dr10 -= mean(dr10,/nan)-mean(fu10,/nan)
        inc10 = dr10-fu10
        
        stop

                                ;115, 162 for June 26
        ;August 6 is bands 94 and 111 (close)
;jul 13 is 276 and 375, also 101 and 111 and 99 , 111
        
;        hy_dr = congrid(hy,n_elements(dr[*,0]),n_elements(dr[0,*]),474,/center)

;first downgrad to 10 m to save conputation cost?
        
;now do NDSI test
        ndsi_c = make_array(474,474,/float,value=nan())
        FOR b1 = 0,473 DO BEGIN
          IF b1 MOD 10 EQ 0 THEN print,'Band ',b1
          FOR b2 = b1,473 DO BEGIN
            IF b1 NE b2 THEN BEGIN
              ndsi = (hy_dr10[*,*,b2]-hy_dr10[*,*,b1])/(hy_dr10[*,*,b2]+hy_dr10[*,*,b1])
              ndsi_c[b1,b2] = correl(inc10,ndsi)
              
            ENDIF 
          ENDFOR
        ENDFOR 

        ndbval = where(~finite(ndsi_c))
        ndsi_c[ndbval] = 0.0
        jj10 = reverse(sort(abs(ndsi_c)))
        FOR jj = 0,10 DO print,array_indices(ndsi_c,jj10[jj]),ndsi_c[jj10[jj]],FORMAT = '(%"Value at [%d, %d] is %f")'

        gband = array_indices(ndsi_c,jj10[0])
        ndsi = (hy_dr10[*,*,gband[1]]-hy_dr10[*,*,gband[0]])/(hy_dr10[*,*,gband[1]]+hy_dr10[*,*,gband[0]])
        incfit = myfit(ndsi,inc10)
        dr_fuh = fu10 + incfit[0] + incfit[1] * ndsi

;save each of these (dr_fu, fu10, dr10, inc10, hy_dr10, ndsi_c,
;jj10, ndsi, incfit, dr_fuh)

        
        stop
        
      ENDIF ELSE BEGIN
        print,'Hyspex and drone do not line up'
        stop
      ENDELSE 
    ENDIF 

    IF keyword_set(lidar) THEN BEGIN
;new way

                                ;extract box for DSM and DEM

      dem_xloc = long(dr_c1[0]-dem_c1[0])
      dem_yloc = long(dem_c1[1]-dr_c1[1])
      dsm_xloc = long(dr_c1[0]-dsm_c1[0])
      dsm_yloc = long(dsm_c1[1]-dr_c1[1])
      
      li_span_x = long(n_elements(dr[*,0])*dr_sc/1.0)
      li_span_y = long(n_elements(dr[0,*])*dr_sc/1.0)

      dem_ext = make_array(li_span_x,li_span_y,/float,value=nan())
      dsm_ext = make_array(li_span_x,li_span_y,/float,value=nan())

      dem_st = 4l * dem_yloc * 30000l
      dsm_st = 4l * dsm_yloc * 30000l
      
      FOR li_ln = 0l,li_span_y-1l DO BEGIN
        dum = fltarr(li_span_x)
        point_lun,dem_fl,dem_st+(4l * li_ln * 30000l)+(dem_xloc*4l)
        readu,dem_fl,dum
        dem_ext[*,li_ln] = dum
        dum = fltarr(li_span_x)
        point_lun,dsm_fl,dsm_st+(4l * li_ln * 30000l)+(dsm_xloc*4l)
        readu,dsm_fl,dum
        dsm_ext[*,li_ln] = dum
      ENDFOR 

      dem_dr = congrid(dem_ext,n_elements(dr[*,0]),n_elements(dr[0,*]),/center)
      dsm_dr = congrid(dsm_ext,n_elements(dr[*,0]),n_elements(dr[0,*]),/center)
      can_dr = congrid(dsm_ext-dem_ext,n_elements(dr[*,0]),n_elements(dr[0,*]),/center)
      
      dem_dr10=congrid(dem_dr,n_elements(dr[*,0])/23,n_elements(dr[0,*])/23,474,/center,/interp)
      dsm_Dr10=congrid(dsm_dr,n_elements(dr[*,0])/23,n_elements(dr[0,*])/23,474,/center,/interp)
      can_dr10=congrid(can_dr,n_elements(dr[*,0])/23,n_elements(dr[0,*])/23,474,/center,/interp)
      stop

;save can_dr, dem_dr, dsm_dr 
      
;      li_loc_x = fix((dr_c1[0]-li_ux)/1.0)
;      li_loc_y = fix((dr_c1[1]-li_uy)/1.0)
;      li_span_x = round(li_loc_x + (n_elements(dr[*,0])*dr_sc/1.0))
;      li_span_y = round(li_loc_y + (n_elements(dr[0,*])*dr_sc/1.0))
;      li_ext = can[li_loc_x:li_span_x,li_loc_y:li_span_y]

      ; test = congrid(can[li_loc_x:li_loc_x+731,li_loc_y+400:li_loc_y+400+530],1575,1258,/center)
    ENDIF 
stop
    
    ;goals: get stats, ,lidar will help?

    ;ecostress overlap days are: 192 (11), 
      
  ENDFOR 

  IF keyword_set(lidar) THEN BEGIN
    free_lun,dsm_fl
    free_lun,dem_fl
  ENDIF 
  

END

PRO ch_dronemodel
  andir = '/air/adesai/cheesehead/lstfusion/analysis/'

  
;restore drone files
  restore,andir+'hy_aug_0712_1630.sav'
  restore,andir+'hy_jul_0712_1630.sav'
  restore,andir+'hy_jun_0712_1630.sav'
;restore lidar files
   restore,andir+'lidr_0712_1630.sav'                              
;read in wwavelenth
   wv = ch_hyspex_wv()


;figure out best ndsi vals across all
;           FOR jj = 0,10 DO
;           print,array_indices(ndsi_c,jj10[jj]),ndsi_c[jj10[jj]],FORMAT = '(%"Value at [%d, %d] is %f")'


;   FOR jj = 0,100 DO print,wv[long((array_indices(ndsi_c_jun,jj10_jun[jj]))[0])],wv[long((array_indices(ndsi_c_jun,jj10_jun[jj]))[1])],ndsi_c_jun[jj10_jun[jj]],ndsi_c_jul[jj10_jun[jj]],ndsi_c_aug[jj10_jun[jj]],format='(5f10.3)'
;   FOR jj = 0,100 DO print,wv[long((array_indices(ndsi_c_jul,jj10_jul[jj]))[0])],wv[long((array_indices(ndsi_c_jul,jj10_jul[jj]))[1])],ndsi_c_jun[jj10_jul[jj]],ndsi_c_jul[jj10_jul[jj]],ndsi_c_aug[jj10_jul[jj]],format='(5f10.3)'
;   FOR jj = 0,100 DO print,wv[long((array_indices(ndsi_c_aug,jj10_aug[jj]))[0])],wv[long((array_indices(ndsi_c_aug,jj10_aug[jj]))[1])],ndsi_c_jun[jj10_aug[jj]],ndsi_c_jul[jj10_aug[jj]],ndsi_c_aug[jj10_aug[jj]],format='(5f10.3)'

;select 1470.466, 1982.657 (281, 375)
;select 709.077 and 760.173 (95, 111)
;select 504.697 and 651.595  (31, 77)

   hy_dr10_jul2 = make_array(68,54,474,/float,value=nan())
   hy_dr10_jul2[*,2:*,*] = hy_dr10_jul
   
   hy_dr10 = merge_array(hy_dr10_jul2,hy_dr10_aug,hy_dr10_jun)
   hy_b3 = (hy_dr10[*,*,375]-hy_dr10[*,*,281])/(hy_dr10[*,*,375]+hy_dr10[*,*,281]) ;SWIR
   hy_b2 = (hy_dr10[*,*,95]-hy_dr10[*,*,111])/(hy_dr10[*,*,95]+hy_dr10[*,*,111])  ;red/green contrast
   hy_b1 = (hy_dr10[*,*,77]-hy_dr10[*,*,31])/(hy_dr10[*,*,77]+hy_dr10[*,*,31])    ;red edge

                                ;fu10, can_dr10, hy_B1, hy_b2, hy_b3

   imx = fltarr(5,68*54)
   imx[0,*] = fu10[*]
   imx[1,*] = can_dr10[*]
   imx[2,*] = hy_b1[*]
   imx[3,*] = hy_b2[*]
   imx[4,*] = hy_b3[*]
   imy = fltarr(1,68*54)
   imy[0,*] = dr10[*]
   imy2 = imy
   imy2[0,*] = inc10[*]
   nn = where(finite(imx[0,*]) AND finite(imx[1,*]) AND finite(imx[2,*]) AND finite(imx[3,*]) AND finite(imx[4,*]) AND finite(imy))
   nn2 = where(finite(imx[1,*]) AND finite(imx[2,*]) AND finite(imx[3,*]) AND finite(imx[4,*]) AND finite(imy2))
   
   t = regress(imx[*,nn],imy[nn],correlation=cor,yfit=yf,const=ct,sigma=sma)
   t2 = regress(imx[2:*,nn2],imy2[nn2],correlation=cor,yfit=yf,const=ct2,sigma=sma)

   ;test combos
   newdr = t[0]*fu10+t[1]*can_dr10+t[2]*hy_b1+t[3]*hy_b2+t[4]*hy_b3+ct
   newdr2 =  t2[0]*hy_b1+t2[1]*hy_b2+t2[2]*hy_b3+ct2+fu10+(mean(dr10,/nan)-mean(fu10,/nan))
   
inctest = t[0]*can_dr10+t[1]*hy_b1+t[2]*hy_b2+t[3]*hy_b3+ct
   
                                ;Aug is usu better than Jul is better
                                ;than Jun (marginally)
stop
   ;take july, merge in aug, then jun

   save,newdr2,hy_b1,hy_b2,hy_b3,fu10,dr10,ndsi_c_jul,ndsi_c_jun,ndsi_c_aug,filename='/air/adesai/cheesehead/lstfusion/analysis/dronemodel.sav'
   
stop
   
;build model
  
END


PRO ch_uwkacomp,version=version
;read IR temp from UWKA

                                ;HOUR
                                ;try TDomeB
                                ;or irbc

                                ;rstb2 (IR radiometric temp) in C

   restore,'/air/adesai/cheesehead/lstfusion/coords.sav'
   uwka_dr = '/bog/incoming/cheesehead/UWKA/'
   files = file_search([uwka_dr+'2019*a.c1.nc',uwka_dr+'2019*b.c1.nc'],count=nf)
   files = files[sort(files)]

   lats = 0.
   lons = 0.
   irs = 0.
   fts = 0.
   gts = 0.
   doys = 0.
   
   FOR i = 0,nf-1 DO BEGIN
     print,'On file ',files[i]
     ncdf_list,files[i],vname=vv,/variables,/quiet
     nc = ncdf_open(files[i])
     dy = strmid((reverse(strsplit(files[i],'/',/extract)))[0],0,8)
     doy = dy_to_jd(dy)
     ncdf_varget,nc,'HOUR',h
     ncdf_varget,nc,'rstb2',ir_temp
     ncdf_varget,nc,'TDomeB',dome_t
     ncdf_varget,nc,'irbc',lw
     IF in(vv,'AVlat') THEN ncdf_Varget,nc,'AVlat',lat
     IF in(vv,'AVlon') THEN ncdf_Varget,nc,'AVlon',lon
     IF in(vv,'avlat') THEN ncdf_Varget,nc,'avlat',lat
     IF in(vv,'avlon') THEN ncdf_Varget,nc,'avlon',lon
     ncdf_varget,nc,'ralt3',z
     ncdf_close,nc

     gv = where((lat LE 46.1 AND lat GE 45.8) AND (lon GE -90.5 AND lon LE -90.0) AND (z GE 50 AND z LE 150),ngv)
     h = h[gv]
     ir_temp = ir_temp[gv]
;     dome_t = dome_t[gv]
;     lw = lw[gv]
     lat = lat[gv]
     lon = lon[gv]
     ir_temp += 273.15

     f_temp = ch_readfusion(lat,lon,sday=doy,eday=doy,doys=replicate(doy,n_elements(lat)),hrs=h,goes=g_temp,/transect,version=version)

     doys = [doys,replicate(doy,n_elements(lat))]
     lats = [lats,lat]
     lons = [lons,lon]
     irs = [irs,ir_temp]
     fts = [fts,f_temp]
     gts = [gts,g_temp]
     
   ENDFOR   

   lats = lats[1:*]
   lons = lons[1:*]
   irs = irs[1:*]
   fts = fts[1:*]
   gts = gts[1:*]
   doys = doys[1:*]
   save,filename='/air/adesai/cheesehead/lstfusion/uwkacomp.sav',irs,fts,gts,doys,lats,lons
stop
   
END



PRO ch_isfscomp2
;read LW_OUT_1_1_1
  rdir = '/bog/incoming/CHEESEHEAD/Ameriflux/'
  dr = '/air/adesai/cheesehead/lstfusion/'
  fname = 'isfs_temperatures.nc'
  version='v20200803'
  nc = ncdf_open(dr+fname)
  ncdf_Varget,nc,'latitude',lat
  ncdf_varget,nc,'longitude',lon
  ncdf_close,nc
  fuse = ch_readfusion(lat,lon,goes=goes,version=version)
  restore,'/air/adesai/cheesehead/lstfusion/emis/aster_em.sav'
  
  st = ['b','c','d','e','g','h','i','j','k','L','m','n','p','q','r','s','t']
  isfs_type = ['ENF','DBF','WET','LAK','ENF','ENF','DBF','ENF','DBF','DBF','DBF','DBF','DBF','DBF','WET','DBF','ENF','LAK']
  sites = strsplit('nw1  nw2  nw3    nw4  ne1  ne2  ne3  ne4  sw1  sw2  sw3  sw4  se2  se3  se4  se5  se6',' ',/extract)
;all sites start 201906010000 LT
  ;end at 201911010000

  tsurf = make_array(17,7344,/float,value=nan())
  tsurf_em = tsurf
  
  FOR i = 0,n_elements(st)-1 DO BEGIN
    fname = rdir+'US-PF'+st[i]+'_HH_201906010000_201911010000.csv'
    print,'Reading ',fname
    dta = read_ascii(fname,delim=',',data_start=1,header=h)
    dta=dta.(0)
    dta[where(Dta LE -9999)]=nan()
    h = strsplit(h,',',/extract)

;18 is LW
                                ;use 1 for emissivity
    tsurf[i,*] = (dta[18,*]/((5.67e-8)*1.0))^0.25

    towerem_x = (closest(em_lon,lon[i]))[0]
    towerem_y = (closest(em_lat,lat[i]))[0]   
    tower_emis = mean(emis[towerem_x-3:towerem_x+3,towerem_y-3:towerem_y+3])
    tsurf_em[i,*] = ( (1/(tower_emis * 5.67e-8)) * (dta[18,*] - ((1-tower_emis)*dta[17,*])))^0.25
;use ASTER GED for emissivity

    
  ENDFOR 


  tsurf_hr = average_cols(tsurf,2,/nan)
  tsurf_hr_em = average_cols(tsurf_em,2,/nan)
  ftemp = shift(fuse,0,-6)
  gtemp = shift(Goes,0,-6)
  stop
  save,filename='/air/adesai/cheesehead/lstfusion/ts.sav',ftemp,gtemp,tsurf_hr,sites,isfs_type,tsurf_hr_em
END


PRO ch_lstfusion2,version=version,keepslope=keepslope,recreate=recreate
;new approach

;one pixel at a time?

                                ;read goesfuse_lst (only 30x30)

                                ;read all 600x600 ecostress (~50 MB)
  indir = '/air/adesai/cheesehead/lstfusion/inputs'
  outdir = '/air/adesai/cheesehead/lstfusion/'
  IF n_elements(version) EQ 0 THEN version = 'v20200803'
  writedir = outdir + 'outputs/'+version+'/'
  create_dir,writedir
  restore,outdir+'goes_lst.sav'
  restore,outdir+'goesfuse_lst.sav'

  IF file_test(writedir+'transform.sav',/read) AND ~keyword_set(recreate) THEN BEGIN
    print,'Restoring ',writedir+'transform.sav'
    restore,writedir+'transform.sav'
    GOTO,runout
  ENDIF 
  
  restore,outdir+'ecostress.sav'

;remove bad images
;  ecostress = screen_arr(Ecostress,271,314)
  neco = n_elements(eco_doy)
  
  ecoblur = make_array(30,30,neco,/float,value=nan())
  FOR i = 0,neco-1 DO BEGIN
    FOR y = 0,29 DO BEGIN
      locy1 = y*20
      locy2 = locy1 + 19
      FOR x = 0,29 DO BEGIN
        locx1 = x*20
        locx2 = locx1 + 19
        ecoblur[x,y,i] = mean(screen_arr(ecostress[locx1:locx2,locy1:locy2,i],271,314),/nan)
      ENDFOR
    ENDFOR
  ENDFOR 
  ecocorr = make_Array(neco,/float,value=nan())
  ecop = ecocorr
  FOR i = 0,neco-1 DO ecocorr[i] = correl(reform(ecoblur[*,*,i]),reform(lsts_g[fix(eco_doy[i])-152,fix(eco_hr[i]),*,*]))
  FOR i = 0,neco-1 DO ecop[i] = corr_ttest(corr=abs(ecocorr[i]),n=900)

  ecokeep = where(ecop lt 0.001 and ecocorr gt 0.3)

  ecostress = ecostress[*,*,ecokeep]
  eco_doy = eco_Doy[ecokeep]
  eco_hr = eco_Hr[ecokeep]
  neco = n_elements(eco_doy)

  lst_out = make_Array(600,600,neco,/float,value=nan())
  
  FOR i = 0,neco-1 DO lst_out[*,*,i] = smooth(congrid(reform(lsts_g[fix(eco_doy[i])-152,fix(eco_hr[i]),*,*]),600,600,/center),[20,20],/edge_mir)

  ecomean = total(ecostress,3,/nan)/total(finite(ecostress),3)

  FOR i = 0,neco-1 DO lst_out[*,*,i] -= ecomean
  FOR i = 0,neco-1 DO ecostress[*,*,i] -= ecomean
  
  slopes = make_array(600,600,/float,value=nan())
  intercepts = make_array(600,600,/float,value=nan())
  corr = make_array(600,600,/float,value=nan())
  pvals = corr
  ecodev = stddev(lst_out-ecostress,/nan)*3
  FOR y = 0,599 DO BEGIN
    IF y MOD 20 EQ 0 THEN print,'Processing row  ',y
    FOR x = 0,599 DO BEGIN
      gv = where(finite(ecostress[x,y,*]) AND (abs(ecostress[x,y,*]-lst_out[x,y,*]) LT ecodev) AND ((ecostress[x,y,*]+ecomean[x,y]) GE 271) AND ((ecostress[x,y,*]+ecomean[x,y]) LE 322),ngv)
      IF (float(ngv)/neco) GT 0.5 THEN BEGIN
        cor = correl(lst_out[x,y,gv],ecostress[x,y,gv])
        corr[x,y] = cor
        pval = corr_ttest(corr=cor,n=ngv)
        pvals[x,y] = pval
        IF pval LT 0.01 AND ngv GT 10 THEN BEGIN
          A = nan()
          B = nan()
          fitexy,lst_out[x,y,gv],ecostress[x,y,gv],A,B,x_sigma=1.,y_sigma=1.
          slopes[x,y] = B
          intercepts[x,y] = A
        ENDIF 
      ENDIF 
    ENDFOR
  ENDFOR 

;gap fill (smooth 1 by 1, solve intercept)

  slp = slopes
  inter = intercepts
  bsl = where(slopes LT 0.9 OR slopes GT 1.5)
  slp[bsl]=nan()
  inter[bsl] = nan()
  mfit = myfit(slopes,intercepts)
  WHILE n_elements(where(~finite(slp))) GT 1 DO slp = merge_array(slp,smooth(slp,[3,3],/nan,/edge_mirror))
  interpred = mfit[0] + mfit[1]*slp
  inter = merge_arraY(inter,interpred)

  save,filename=writedir+'transform.sav',slopes,intercepts,slp,inter,corr,pvals,lst_out,ecostress,ecocorr,ecop,ecokeep,eco_doy,eco_hr,neco,ecomean

  
  stop

  
                                ;'fuse_YYYY_DOY_HH.dat'
                                ;lsts_g = [153,24,*,*]
  runout:

  IF keyword_set(keepslope) THEN BEGIN
    mfit = myfit(slopes,intercepts)
    WHILE n_elements(where(~finite(slopes))) GT 1 DO slopes = merge_array(slopes,smooth(slopes,[3,3],/nan,/edge_mirror))
    interpred = mfit[0] + mfit[1]*slopes
    intercepts = merge_arraY(intercepts,interpred)
   ENDIF 
    
  FOR doy = 152,304 DO BEGIN
    print,'Writing day ',doy,' of 304'
    FOR hr = 0,23 DO BEGIN
      thelst = smooth(congrid(reform(lsts_g[doy-152,hr,*,*]),600,600,/center),[20,20],/edge_mir)
      IF keyword_set(keepslope) THEN newlst = (thelst-ecomean)*slopes+intercepts+ecomean ELSE newlst = (thelst-ecomean)*slp+inter+ecomean
      img = fix(newlst*10)      
      doystr = string(doy,format='(i3.3)')
      hrstr = string(Hr,format='(i2.2)')
      openw,fl,writedir+'fuse_2019_'+doystr+'_'+hrstr+'.dat',/get_lun
      writeu,fl,img
      free_lun,fl
    ENDFOR
  ENDFOR 

;output for all hours; call it good

                                ;re-run utm and drone comparison

                                ;redo domain for gridding - UTM

                                ;NLDAS for outer, fusion for inner

  ;later: add lidar
  
;for each 600x600 - 
;correlate (fitexy)  x = goes y = ecostress,

                                ;save m and b

                                ;then for each file - write out goes_lst * ms + bs
  ;perhaps do some spatial filtering?

  

END

PRO ch_fusiondev
  fdev = make_array(24,153,/float,value=nan())
  gdev = fdev
  outdir = '/air/adesai/cheesehead/lstfusion/'
  version = 'v20200803'
  writedir = outdir + 'outputs/'+version+'/'
  restore,outdir+'goesfuse_lst.sav'
  FOR doy = 152,304 DO BEGIN
    print,'Read day ',doy,' of 304'
    FOR hr = 0,23 DO BEGIN
      doystr = string(doy,format='(i3.3)')
      hrstr = string(Hr,format='(i2.2)')
      openr,fl,writedir+'fuse_2019_'+doystr+'_'+hrstr+'.dat',/get_lun
      img = intarr(600,600)
      readu,fl,img
      free_lun,fl
      img = float(img)/10.
      fdev[hr,doy-152]=stddev(img,/nan)
      gdev[hr,doy-152]=stddev(lsts_g[doy-152,hr,*,*])
    ENDFOR
  ENDFOR 

stop
  
END


;read all

PRO ch_readalleco
  g = file_search('/air/adesai/cheesehead/lstfusion/ecostress/ecostress_*.dat',count=ng)
  ecostress = make_array(600,600,ng,/float,value=nan())
  eco_doy = fltarr(ng)
  eco_hr = fltarr(ng)
  FOR i = 0,ng-1 DO BEGIN
    print,'reading ',g[i]
    openr,fl,g[i],/get_lun
    ec = intarr(600,600)
    readu,fl,ec
    freE_lun,fl
    ec = float(ec)
    ec[where(Ec LE -999)]=nan()
    ec/=10.
    ecostress[*,*,i]=ec
    eco_doy[i] = float(strmid((reverse(strsplit(g[i],'/',/extract)))[0],15,3))
    eco_hr[i] = float(strmid((reverse(strsplit(g[i],'/',/extract)))[0],19,2))
  ENDFOR 

  save,ecostress,eco_doy,eco_hr,filename='/air/adesai/cheesehead/lstfusion/ecostress.sav'
  eco_tm = eco_doy + eco_hr/24.
  stop
END


PRO ch_lstout,version=version
;sreenath output
;-90.538400,45.765199, 800
;-90.022154,45.765199,
;-90.021371,46.125030,
;-90.539185,46.125030.

;read NLDAS and Fusion

  outdir = '/air/adesai/cheesehead/lstfusion/'
  IF n_elements(version) EQ 0 THEN version = 'v20200803'
  fusedir = outdir + 'outputs/'+version+'/'
  writedir = outdir + 'geotiff/'+version+'/'
  create_dir,writedir

  inlat = numgen(46.1,45.8,n=600)
  inlon = numgen(-90.5,-90,n=600)

  outlat = numgen(46.125030,45.765199,n=800)
  outlon = numgen(-90.539185,-90.021371,n=800)

  glat = where(outlat LE inlat[0] AND outlat GE inlat[599],nglat)
  glon = where(outlon GE inlon[0] AND outlon LE inlon[599],nglon)

  outlat_nl = numgen(46.125030,45.765199,n=9)
  outlon_nl = numgen(-90.539185,-90.021371,n=9)

  geo = {ModelTiepointTag : [0,0, 0, outlon[0],outlat[0], 0],  ModelPixelScaleTag  : [outlon[1]-outlon[0],outlat[0]-outlat[1],1.0], $
         GTModelTypeGeoKey : 2, GTRasterTypeGeoKey : 2, GEOGRAPHICTYPEGEOKEY : 4326, GEOCITATIONGEOKEY : 'GCS_WGS_84'}    
  

  FOR doy = 152,304 DO BEGIN
    print,'Day ',doy
    doystr = string(doy,format='(i3.3)')
    ch_nldas,doys=doy,/nosave,outlat=outlat_nl,outlon=outlon_nl,output=nldas
    FOR hr = 0,23 DO BEGIN
      print,'  Hour ',hr
      hrstr = string(hr,format='(i2.2)')
      img = intarr(600,600)
      openr,fl,fusedir+'fuse_2019_'+doystr+'_'+hrstr+'.dat',/get_lun
      readu,fl,img
      free_lun,fl
      img = img / 10.      
      nld = congrid(total(reform(nldas[hr,*,*,*]),3)/3.0,800,800,/interp)
      nld = nld + mean(img-nld,/nan)
      outarr = nld
      FOR lat = 0,nglat-1 DO BEGIN
        FOR lon = 0,nglon-1 DO BEGIN
          nld[glon[lon],glat[lat]] = img[closest(inlon,outlon[glon[lon]]),closest(inlat,outlat[glat[lat]])]          
        ENDFOR
      ENDFOR
      nld2 = smooth(nld,[15,15],/edge_mir)
      nld2[glon[0]:glon[nglon-1],glat[0]:glat[nglat-1]] = nld[glon[0]:glon[nglon-1],glat[0]:glat[nglat-1]]
      write_Tiff,writedir+'fuse_2019_'+doystr+'_'+hrstr+'.tif',nld2,/float,geotiff=geo,compress=1
    ENDFOR 
  ENDFOR 

  ;geotiff

END


PRO ch_asteremis
  a1 = '/air/adesai/cheesehead/lstfusion/emis/AG100.v003.47.-091.0001.h5'
  a2 = '/air/adesai/cheesehead/lstfusion/emis/AG100.v003.46.-091.0001.h5'

  em1 = h5_getdata(a1,'/Emissivity/Mean')
  em2 = h5_getdata(a2,'/Emissivity/Mean')

  em13 = [[em1[*,*,3]],[em2[*,*,3]]]
  em14 = [[em1[*,*,4]],[em2[*,*,4]]]

                                ;46.1 to 45.8 = 900 to 1200
                                ;-90.5 to 90 = 500 to 1000
  em13 = em13[500:999,900:1199]/1000.
  em14 = em14[500:999,900:1199]/1000.
  
  em_lat = numgen(46.1,45.8,n=300)
  em_lon = numgen(-90.5,-90,n=500)

  e10 = 0.6820 + 0.2578 * em13 + 0.0584 * em14
  e11 = -.5415 + 1.4305 * em13 + 0.1092 * em14
  emis = (e10+e11)/2.0

  save,filename='/air/adesai/cheesehead/lstfusion/emis/aster_em.sav',em_lat,em_lon,emis

END


PRO ch_lst_figures
;study area with grids, images
  dr = '/air/adesai/cheesehead/lstfusion/'
  version = 'v20200803/'


  restore,'/air/adesai/cheesehead/lstfusion/ts.sav'
  restore,dr+'goesfuse_lst.sav'
  restore,dr+'goes_lst.sav'
  restore,dr+'goesfuse_lst_stats.sav'
  restore,dr+'nldas_lst.sav'
  restore,dr+'outputs/'+version+'transform.sav'
  restore,dr+'analysis/fdev.sav'
  glat = numgen(46.1,45.8,n=600)
  glon = numgen(-90.5,-90,n=600)
  avsft_m = total(avsft_i,5,/nan)/total(finite(avsft_i),5)
  avsft_s = stddev(avsft_i,dim=5)
                                ;sday=152, eday=304

;idea variance of 
  !p.multi = [0,2,2]
  loadct,49
  tvplot,fdev-gdev,zrange=[-0.2,1.3],xrange=[0,24],yrange=[153,204],/order,ytitle='Day of Year',xtitle='Hour (UTC)'
  
;FIGURE A (read images, plot a,b,c. skip drone nd hyspex)
;ecostress image, goes image, NLDAS, + map
;day 219 0Z
  !p.multi = [0,2,2]
  plot,[0,0],[0,0],/nodata,title='a) 7 Aug 2019 0Z'
  loadct,70
  tvplot,avsft_m[67,0,*,*],zrange=[298,291],title='b)',xrange=[min(glon),max(glon)],yrange=[min(glat),max(glat)],charsize=2
  tvplot,lsts_i[67,0,*,*],zrange=[298,291],title='c)',xrange=[min(glon),max(glon)],yrange=[min(glat),max(glat)],charsize=2
  tvplot,ecostress[*,*,15]+ecomean,xrange=[min(glon),max(glon)],yrange=[min(glat),max(glat)],charsize=2,title='d)',zrange=[298,291]
  
  stop

  golat = numgen(46.1,45.8,n=30)
  golon = numgen(-90.5,-90,n=30)

                                ;FIGURE 4
  
  ;SE3 wetland - 13
  lat_a = 45.927149999999997
  lon_a = -90.247500000000002
  temp_a = ch_readfusion(lat_a,lon_a,version='v20200803')
  glat1 = closest(golat,lat_a)
  glon1 = closest(golon,lon_a)
  pldoy = numgen(152.,304.,numvals=153*24l)
  nl_y = (transpose(reform(avsft_m[*,*,glon1,glat1])))[*]
  nl_ye = (transpose(reform(avsft_s[*,*,glon1,glat1])))[*]
  tower_y = shift(reform(tsurf_hr_em[13,*]),6)
  
   ;SE4 DBF
  lat_b = 45.924483299999999
  lon_b = -90.247450000000001
  
  temp_b = ch_readfusion(lat_b,lon_b,version='v20200803')
  glat2 = closest(golat,lat_b)
  glon2 = closest(golon,lon_b)
  nl_yb = (transpose(reform(avsft_m[*,*,glon2,glat2])))[*]
  nl_yeb = (transpose(reform(avsft_s[*,*,glon2,glat2])))[*]
  tower_yb = shift(reform(tsurf_hr_em[14,*]),6)
  
  ;SE5 ENF
  lat_c = 45.938083300000002
  lon_c = -90.238183300000003
  temp_c = ch_readfusion(lat_c,lon_c,version='v20200803')
  glat3 = closest(golat,lat_c)
  glon3 = closest(golon,lon_c)
  nl_yc = (transpose(reform(avsft_m[*,*,glon3,glat3])))[*]
  nl_yec = (transpose(reform(avsft_s[*,*,glon3,glat3])))[*]
  tower_yc = shift(reform(tsurf_hr_em[15,*]),6)
  
   ;SE6 LAK
  lat_d = 45.919733299999997
  lon_d = -90.228833300000005
  temp_d = ch_readfusion(lat_d,lon_d,version='v20200803')
  glat4 = closest(golat,lat_d)
  glon4 = closest(golon,lon_d)
  nl_yd = (transpose(reform(avsft_m[*,*,glon4,glat4])))[*]
  nl_yed = (transpose(reform(avsft_s[*,*,glon4,glat4])))[*]
  tower_yd = shift(reform(tsurf_hr_em[16,*]),6)

  !p.multi = [0,1,4]
  red = fsc_color('red',201)
  blue = fsc_color('blue',202)
  orange = fsc_color('orange',203)
  
;se3
  plot,[0,0],[0,0],/nodata,xrange=[195,205],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  errorbar,pldoy,nl_y,nl_ye,nl_ye,fillcolor=200,line_fill=0,linecolor=0,/over
  oplot,pldoy,nl_y,thick=2
  oplot,pldoy,transpose(reform(lsts_i[*,*,glon1,glat1])),thick=2,color=blue,psym=-1
  oplot,pldoy,temp_a,thick=2,color=red
  oplot,pldoy,tower_y,thick=3,color=orange
  ;; plot,[0,0],[0,0],/nodata,xrange=[270,275],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  ;; errorbar,pldoy,nl_y,nl_ye,nl_ye,fillcolor=200,line_fill=0,linecolor=0,/over
  ;; oplot,pldoy,nl_y,thick=2
  ;; oplot,pldoy,transpose(reform(lsts_i[*,*,glon1,glat1])),thick=2,color=blue,psym=-1
  ;; oplot,pldoy,temp_a,thick=2,color=red
  ;; oplot,pldoy,tower_y,thick=2,color=orange
  
;se4
  plot,[0,0],[0,0],/nodata,xrange=[195,205],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  errorbar,pldoy,nl_yb,nl_yeb,nl_yeb,fillcolor=200,line_fill=0,linecolor=0,/over
  oplot,pldoy,nl_yb,thick=2
  oplot,pldoy,transpose(reform(lsts_i[*,*,glon2,glat2])),thick=2,color=blue,psym=-1
  oplot,pldoy,temp_b,thick=2,color=red
  oplot,pldoy,tower_yb,thick=3,color=orange
  ;; plot,[0,0],[0,0],/nodata,xrange=[270,275],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  ;; errorbar,pldoy,nl_yb,nl_yeb,nl_yeb,fillcolor=200,line_fill=0,linecolor=0,/over
  ;; oplot,pldoy,nl_yb,thick=2
  ;; oplot,pldoy,transpose(reform(lsts_i[*,*,glon2,glat2])),thick=2,color=blue,psym=-1
  ;; oplot,pldoy,temp_b,thick=2,color=red
  ;; oplot,pldoy,tower_yb,thick=2,color=orange

;se5
  plot,[0,0],[0,0],/nodata,xrange=[195,205],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  errorbar,pldoy,nl_yc,nl_yec,nl_yec,fillcolor=200,line_fill=0,linecolor=0,/over
  oplot,pldoy,nl_yc,thick=2
  oplot,pldoy,transpose(reform(lsts_i[*,*,glon3,glat3])),thick=2,color=blue,psym=-1
  oplot,pldoy,temp_c,thick=2,color=red
  oplot,pldoy,tower_yc,thick=3,color=orange
  ;; plot,[0,0],[0,0],/nodata,xrange=[270,275],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  ;; errorbar,pldoy,nl_yc,nl_yec,nl_yec,fillcolor=200,line_fill=0,linecolor=0,/over
  ;; oplot,pldoy,nl_yc,thick=2
  ;; oplot,pldoy,transpose(reform(lsts_i[*,*,glon3,glat3])),thick=2,color=blue,psym=-1
  ;; oplot,pldoy,temp_c,thick=2,color=red
  ;; oplot,pldoy,tower_yc,thick=2,color=orange
    
;se6
  plot,[0,0],[0,0],/nodata,xrange=[195,205],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  errorbar,pldoy,nl_yd,nl_yed,nl_yed,fillcolor=200,line_fill=0,linecolor=0,/over
  oplot,pldoy,nl_yd,thick=2
  oplot,pldoy,transpose(reform(lsts_i[*,*,glon4,glat4])),thick=2,color=blue,psym=-1
  oplot,pldoy,temp_d,thick=2,color=red
  oplot,pldoy,tower_yd,thick=3,color=orange
  ;; plot,[0,0],[0,0],/nodata,xrange=[270,275],yrange=[270,310],xtitle='DOY',ytitle='Degrees K',title='a)',charsize=1.25
  ;; errorbar,pldoy,nl_yd,nl_yed,nl_yed,fillcolor=200,line_fill=0,linecolor=0,/over
  ;; oplot,pldoy,nl_yd,thick=2
  ;; oplot,pldoy,transpose(reform(lsts_i[*,*,glon4,glat4])),thick=2,color=blue,psym=-1
  ;; oplot,pldoy,temp_d,thick=2,color=red
  ;; oplot,pldoy,tower_yd,thick=2,color=orange
    

  stop
  
  ;FIGURE B
;GOES-NLDAS scatter plot

  gn_h = hist_2d(nlfix[*],lsts_i[*],bin1=1,bin2=1,min1=250,min2=250,max1=310,max2=310)
  tvplot,alog10(gn_h>1),/order,zrange=[5,0],charsize=2,xrange=[250,310],yrange=[250,310]
  oplot,[0,1000],[0,1000],thick=2
                                ;r=0.905239 (range r^2 70-88% var)
                                ;RMSE=3.58
  ;bias = -0.78 K

  ;original NLDAS had r=0.895 RMSE = 3.68
  stop
  

  ;FIGURE C
                                ;ECOSTRESS-GOES correlation
                                ;see lstfusion2
  ;restore,dr+'ecostress.sav'


  !p.multi=[0,3,3]
  loadct,70
  tvplot,corr,xrange=[min(glon),max(glon)],yrange=[min(glat),max(glat)],charsize=2,zrange=[1,0.7],title='Correlation'
  tvplot,slp,xrange=[min(glon),max(glon)],yrange=[min(glat),max(glat)],charsize=2,zrange=[1.5,0.7],title='Slope'
  tvplot,inter,xrange=[min(glon),max(glon)],yrange=[min(glat),max(glat)],charsize=2,zrange=[3,-3],title='Intercept'
  
;transform has slopes, intercepts, corrected: slp, inter
                                ;corr, pvals, lst_out (ecostress) ecoblur - upscaled lsts_g
  ;ecostress, ecocorr,ecop (keep good images - how many?)
;25 images out of 49 were retained with r from 0.32 to 0.74 (P<0.001)

  ;temporal corre is bet 0.59 ti 0.95, pvalues from 0 to 0.012 (P<0.01)
   

  stop


  ;FIGURE D (figure 4) DO THIS ONE
                                ;FUSION zoom in - diel cycle,
                                ;lake, forest, city
                                ;GOES, FUSE, 5 days in june and sept

  ;pick spots to do diurnals from lsts_g, nldas, and fuse

  

  ;FIGURE E
;Flux tower and UWKA comparisons isfs_comp and uwka_comp
;tower color by hour, symbol by type
;isfs - 

  restore,'/air/adesai/cheesehead/lstfusion/uwkacomp.sav'
  restore,'/air/adesai/cheesehead/lstfusion/ts.sav'

  uwka_hist = hist_2d(fts,irs,bin1=1,bin2=1,min1=280,min2=280,max1=315,max2=315)

  
  !p.multi = [0,2,2]
  tvplot,alog10(uwka_hist>1),/order,zrange=[5,0],charsize=2,xrange=[280,315],yrange=[280,315]

  dbf = fsc_color('green',200)
  enf = fsc_Color('brown',201)
  wet = fsc_color('cyan',202)
  isfs_col = [enf,dbf,wet,wet,enf,enf,dbf,dbf,dbf,dbf,dbf,dbf,dbf,dbf,wet,dbf,enf]
  plot,[0,0],[0,0],/nodata,xrange=[260,310],yrange=[260,310],xtitle='Fusion Degrees K',ytitle='Tower Degrees K',charsize=2
  FOR i = 200,202 DO oplot,ftemp[where(isfs_col EQ i),*],tsurf_hr_em[where(isfs_col EQ i),*],color=i,psym=1,symsize=0.5
  oplot,[200,400],[200,400],thick=2
  stop
  
                                ;uwka r = 0.864 for fusion, 0.856 for
                                ;goes. RMSE is 2.49 K

;by individual days r varies from 0.28 to 0.82, mean 0.55
  
                                ;towers, r = 0.875, RMSE is 4.2 K
                                ;correlation for enf and dbf are both
                                ;0.89 dbf, 0.88 for enf, but 0.85 for
                                ;wetland

  

;uwka - color by IOP

  

  
  ;FIGURE F
;Landsat - read/compare
  restore,'/air/adesai/cheesehead/lstfusion/analysis/landsat.sav'
  !p.multi=[0,2,2]
  tvplot,lsat_rg[0:400,*],xrange=[694300l,694300l+400l*50l],yrange=[5075900l,5075900l+668l*50l],xticks=1,yticks=1,xtickformat='(i0)',ytickformat='(i0)',zrange=[285,293],charsize=1.5,xtitle='Easting',ytitle='Northing'
  tvplot,fuse_lsat[0:400,*],xrange=[694300l,694300l+400l*50l],yrange=[5075900l,5075900l+668l*50l],xticks=1,yticks=1,xtickformat='(i0)',ytickformat='(i0)',zrange=[285,293],charsize=1.5,xtitle='Easting',ytitle='Northing'
  lsat_hist = hist_2d(lsat_Rg,fuse_lsat,bin1=1,bin2=1,min1=285,min2=285,max1=295,max2=295)
  tvplot,alog10(lsat_hist>1),/order,zrange=[5,0],charsize=1.5,xrange=[285,295],yrange=[285,295],xtitle='Landsat',ytitle='Fusion'
  oplot,[0,1000],[0,1000],thick=2
  stop
  
  ;r = 0.52

  ;FIGURE G
;Drone model - image, show spectra for sample areas, and lines

  drone_dr = '/bog/incoming/cheesehead/drone-lst/'
  flights = ['07112019','07112019','07112019',$
             '07122019','07122019','07122019','07122019',$
             '08192019','08192019',$
             '08212019','08212019','08212019','08212019','08212019']
  fnumber = ['1','2','4','1','2','3','4','1','2','1','2','3','4','5']
  dr = read_tiff(drone_dr+flights[6]+'_flight'+fnumber[6]+'.tif',geotiff=gdr)
  dr[where(dr LT 0)]=nan()
  dr+=273.15
  restore,'/air/adesai/cheesehead/lstfusion/analysis/dronemodel.sav'
  wv = ch_hyspex_wv()

                                ;wv is not linear
  newwv = numgen(405,2517,3.)
  ndsi_inter = fltarr(n_elements(newwv),n_elements(newwv))
  FOR i = 0,n_elements(newwv)-1 DO FOR j = 0,n_elements(newwv)-1 DO ndsi_inter[i,j] = ndsi_c_aug[closest(wv,newwv[i]),closest(wv,newwv[j])]


;NDSI plot
  loadct,70
  tvplot,ndsi_inter,xrange=[min(newwv),max(newwv)],yrange=[min(newwv),max(newwv)],zrange=[0.5,-0.5],xtitle='nm',ytitle='nm',/order
  stop 
  
  dr_sc = gdr.modelpixelscaletag[0]
  dr_c1 = gdr.modeltiepointtag[3:4]
  dr_east = dr_c1[0]+dr_sc*n_elements(dr[*,0])
  dr_north = dr_c1[1]-dr_sc*n_elements(dr[0,*])
  
  !p.multi=[0,2,2]
                                ;dr dr10 dr_fu dr2
  tvplot,dr,zrange=[295,320],xrange=[dr_c1[0],dr_east],yrange=[dr_north,dr_C1[1]],xtickformat='(i0)',xticks=1,ytickformat='(i0)',yticks=1
  tvplot,dr10,zrange=[295,320],xrange=[dr_c1[0],dr_east],yrange=[dr_north,dr_C1[1]],xtickformat='(i0)',xticks=1,ytickformat='(i0)',yticks=1
  tvplot,fu10,zrange=[295,320],xrange=[dr_c1[0],dr_east],yrange=[dr_north,dr_C1[1]],xtickformat='(i0)',xticks=1,ytickformat='(i0)',yticks=1
  tvplot,newdr2,zrange=[295,320],xrange=[dr_c1[0],dr_east],yrange=[dr_north,dr_C1[1]],xtickformat='(i0)',xticks=1,ytickformat='(i0)',yticks=1

  ;before downscale r = 0.37, after downscale r = 0.58
  
  stop

;FIGURE H
  ;power spectrum UWKA vs FUSE vs GOES
  
; pirs2 =
; fft_powerspectrum(smooth(irs,10),freq=ff,/tukey,width=0.01,sign=signif)
                                ;average by doy

; freqDomainImage = fft(k,-1)
  ;take nlfix vs ecostress vs read lstfusion
;power = SHIFT(ALOG(ABS(freqDomainImage)), xsize/2, ysize/2)


                                ;600x600x25
                                ;4 and 15 are nice
                                ;4 = 167, 20; 15 = 216, 0

  fu_comp = '/air/adesai/cheesehead/lstfusion/outputs/'+version+'fuse_2019_167_20.dat'
  fu_compd = intarr(600,600)
  openu,fl,fu_comp,/get_lun
  readu,fl,fu_compd
  free_lun,fl
  fu_compd/=10.

  ec = reform(ecostress[*,*,4])
  ec = zapbadval(ec)
  ec+=273.15

  nl_comp = reform(nlfix[15,20,*,*])
  
  power_fu = shift(alog(abs(fft(fu_compd,-1))), 300, 300)
  power_ec = shift(alog(abs(fft(ec,-1))),300,300)
  power_nl = shift(alog(abs(fft(nl_comp,-1))),15,15)

  loadct,33
  !p.multi =[0,3,3]
  tvplot,power_nl,zrange=[-14,5],xrange=[-30,30],yrange=[-30,30]
  tvplot,power_ec,zrange=[-14,5],xrange=[-30,30],yrange=[-30,30]
  tvplot,power_fu,zrange=[-14,5],xrange=[-30,30],yrange=[-30,30]


  stop
  
  
END

PRO ch_landsat
  landsatdr = '/air/incoming/cheesehead/landsatLST/'
  drs= ['LC08_L1TP_025028_20190615_20190620_01_T1',$
        'LC08_L1TP_026028_20190708_20190719_01_T1',$
        'LC08_L1TP_026028_20190724_20190801_01_T1',$
        'LC08_L1TP_026028_20190809_20190820_01_T1',$
        'LC08_L1TP_026028_20190926_20191017_01_T1']

  yyyymmdd = ['20190615','20190708','20190724','20190809','20190926']
  doy = dy_to_jd(yyyymmdd)
  
  fnames = landsatdr + drs + '/' + drs + '_SW_LST.tif'
  hr = 16
  vdir = 'v20200803/'
  version = 'v20200803'
    
  i = 4  ;2 is good so far, 3 ok, 4 is good (r=0.51)
  lsat = read_tiff(fnames[i],geotiff=gls)
  sc = gls.modelpixelscaletag[0]
  IF i EQ  0 THEN BEGIN
    tp0 = gls.modeltiepointtag[3:4]
    tpll = utm_to_ll(tp0[0],tp0[1],'16N',Ref='WGS84')
    tp = ll_to_utm(tpll[0],tpll[1])
  ENDIF ELSE BEGIN 
;30 m
    tp = gls.modeltiepointtag[3:4]
  ENDELSE 


  
;0
                                ;47.06447, -91.55429

  ;landsat is 10 am (16:46) so 16:30
  ;tp0 = ll_to_utm(47.06447, -91.55429)
  
  x1 = (694300-tp[0])/sc
  x2 = (731600-tp[0])/sc
  y1 = (tp[1]-5109300)/sc
  y2 = (tp[1]-5075900)/sc

  lsat_ex = lsat[x1:x2,y1:y2]
  lsat_rg = congrid(lsat_ex,747,669,/center)
  
  
                                ;extract 747 by 669
  fname = '/air/adesai/cheesehead/lstfusion/outputs/'+vdir+'utm/fuse_2019_'+string(doy[i],format='(i3.3)')+'_'+string(hr,format='(i2.2)')+'.tif'
  IF ~file_test(fname,/read) THEN ch_reproject_write,doy[i],hr,version=version
  IF file_test(fname,/read) THEN fuse = read_tiff(fname,geotiff=gfu)


  
stop

lsat_rg[where(lsat_rg LT min(fuse))]=nan()
lsat_doy = doy[i]
fuse_lsat = fuse
save,filename='/air/adesai/cheesehead/lstfusion/analysis/landsat.sav',lsat_rg,lsat_doy,fuse_lsat

;   694300 to 731600 east
; 5109300 5075900 north

END



;code to read hyspex images

;identify files

FUNCTION ch_hyspex_hdr,fname
  openr,fl,fname,/get_lun
  s=''
  readf,fl,s
  IF s EQ 'ENVI' THEN BEGIN
    WHILE ~eof(fl) DO BEGIN
      s = ''
      readf,fl,s
      IF strmid(s,0,5) EQ 'lines' THEN height = long(strmid(s,8))
      IF strmid(s,0,7) EQ 'samples' THEN width = long(strmid(s,10))
      IF strmid(s,0,5) EQ 'bands' THEN bands = long(strmid(s,8))
      IF strmid(s,0,9) EQ 'data type' THEN datatype = long(strmid(s,12))
      IF strmid(s,0,8) EQ 'map info' THEN BEGIN
        r = strsplit(strmid(s,12),',',/extract)
        east = float(r[3])
        north = float(r[4])
        resolution = float(r[5])
      ENDIF 
    ENDWHILE
    free_lun,fl
    return,[east,north,width,height,bands,resolution,datatype]
  ENDIF ELSE BEGIN
    print,'Non ENVI file ',fname
    return,[-1,-1,0,0,0,0,0]
  ENDELSE 
END

PRO ch_hyspex_find,dr,easting,northing,files=files,coords=coords,locs=locs,output=output,fls=fls
;return list of files that are in points

  fls = file_search(dr+'*.hdr',count=nf)

  IF nf GT 0 THEN BEGIN
    print,'Found ',nf,' files'
    output = make_Array(7,nf,/float,value=nan())
    FOR i = 0,nf-1 DO output[*,i]=ch_hyspex_hdr(fls[i])
    ncoord = n_elements(easting)
    files = strarr(ncoord)
    coords = make_Array(7,ncoord,/float,value=nan())
    locs = make_array(2,ncoord,/long,value=nan())
    e1 = output[0,*]
    e2 = output[0,*] + (output[2,*]*output[5,*])
    n3 = output[1,*]
    n4 = output[1,*] - (output[3,*]*output[5,*])
    FOR i = 0,ncoord-1 DO BEGIN
      ccheck = ((easting[i]-e1) GE 0) AND $
               ((e2-easting[i]) GE 0) AND $
               ((n3-northing[i]) GE 0) AND $
               ((northing[i]-n4) GE 0)
      cloc = where(ccheck,ncl)
      IF ncl GT 0 THEN BEGIN
        print,'Found ',ncl,' matches for coordinate ',i
        print,'  ',fls[cloc]
        cloc = cloc[ncl/2]
        coords[*,i] = output[*,cloc]
        files[i] = fls[cloc]
        locs[0,i] = long((easting[i]-e1[cloc])/coords[5,i])
        locs[1,i] = long((n3[cloc]-northing[i])/coords[5,i])
      ENDIF 
    ENDFOR     
  ENDIF ELSE BEGIN
    print,'No files found matching '+dr+'*.hdr'
  ENDELSE
END

FUNCTION ch_hyspex_rd,fl,box=box,bands=bands
                                ;given file, open,  read, subset
  IF file_test(fl,/read) AND file_test(fl+'.hdr',/read) THEN BEGIN
    hdr = ch_hyspex_hdr(fl+'.hdr')
    IF hdr[0] NE -1 THEN BEGIN
      IF n_elements(bands) EQ 0 THEN bands = 0
      cols = long64(hdr[2])
      rows = long64(hdr[3])
      datatype = long64(hdr[6])
      IF n_elements(box) EQ 0 THEN box = [0,0,cols-1,rows-1]
      ncols = long64((box[2]-box[0])+1)
      nrows = long64((box[3]-box[1])+1)

      
      output = make_array(ncols,nrows,n_elements(bands),/float,value=nan())
      openr,f,fl,/get_lun
      FOR b = long64(0),long64(n_elements(bands)-1) DO BEGIN
        FOR r = long64(0),long64(nrows-1) DO BEGIN
          IF datatype EQ 4 THEN info = fltarr(ncols) ELSE info = intarr(ncols)
          point_lun,f,datatype*((cols*rows*b)+((box[1]+r)*cols)+box[0])
          readu,f,info
          IF datatype EQ 2 THEN info/=10000.
          output[*,r,b] = info
        ENDFOR 
      ENDFOR
      free_lun,f
      return,output
    ENDIF ELSE BEGIN
      print,'Stopping read'
      return,-1
    ENDELSE 
  ENDIF ELSE BEGIN
    print,'Cannot find file and/or header ',fl
    return,-1
  ENDELSE

                                ;if data type = 2, int16, data type = 4 float32
  ;if int16 then refl = refl*10000
    
END

FUNCTION ch_hyspex_wv
  wavelength = [405.7000607383701, 408.89350694882324, 412.08695315927645, 415.2803993697296, 418.4738455801828 $
              , 421.66729179063594, 424.86073800108915, 428.0541842115423, 431.24763042199544, 434.44107663244864 $
              , 437.6345228429018, 440.827969053355, 444.02141526380814, 447.21486147426134, 450.4083076847145 $
              , 453.60175389516763, 456.79520010562084, 459.988646316074, 463.1820925265272, 466.37553873698033 $
              , 469.56898494743353, 472.7624311578867, 475.9558773683399, 479.14932357879303, 482.34276978924623 $
              , 485.5362159996994, 488.7296622101525, 491.9231084206057, 495.1165546310589, 498.310000841512 $
              , 501.5034470519652, 504.69689326241837, 507.8903394728716, 511.0837856833247, 514.2772318937779 $
              , 517.4706781042311, 520.6641243146842, 523.8575705251374, 527.0510167355906, 530.2444629460438 $
              , 533.437909156497, 536.6313553669501, 539.8248015774033, 543.0182477878565, 546.2116939983096 $
              , 549.4051402087628, 552.598586419216, 555.7920326296692, 558.9854788401224, 562.1789250505755 $
              , 565.3723712610285, 568.5658174714818, 571.759263681935, 574.9527098923882, 578.1461561028414 $
              , 581.3396023132944, 584.5330485237477, 587.7264947342007, 590.919940944654, 594.1133871551071 $
              , 597.3068333655604, 600.5002795760136, 603.6937257864666, 606.8871719969198, 610.0806182073729 $
              , 613.2740644178261, 616.4675106282793, 619.6609568387325, 622.8544030491857, 626.0478492596388 $
              , 629.241295470092, 632.4347416805451, 635.6281878909983, 638.8216341014515, 642.0150803119047 $
              , 645.2085265223579, 648.401972732811, 651.5954189432642, 654.7888651537173, 657.9823113641705 $
              , 661.1757575746237, 664.3692037850769, 667.56264999553, 670.7560962059832, 673.9495424164364 $
              , 677.1429886268895, 680.3364348373427, 683.5298810477959, 686.7233272582491, 689.9167734687022 $
              , 693.1102196791554, 696.3036658896086, 699.4971121000617, 702.6905583105149, 705.8840045209681 $
              , 709.0774507314213, 712.2708969418744, 715.4643431523276, 718.6577893627808, 721.8512355732339 $
              , 725.0446817836871, 728.2381279941403, 731.4315742045935, 734.6250204150466, 737.8184666254998 $
              , 741.011912835953, 744.2053590464061, 747.3988052568593, 750.5922514673125, 753.7856976777657 $
              , 756.9791438882188, 760.172590098672, 763.3660363091252, 766.5594825195783, 769.7529287300315 $
              , 772.9463749404847, 776.1398211509379, 779.333267361391, 782.5267135718442, 785.7201597822974 $
              , 788.9136059927505, 792.1070522032037, 795.3004984136569, 798.4939446241101, 801.6873908345632 $
              , 804.8808370450164, 808.0742832554696, 811.2677294659227, 814.4611756763759, 817.6546218868291 $
              , 820.8480680972823, 824.0415143077354, 827.2349605181886, 830.4284067286418, 833.6218529390949 $
              , 836.8152991495481, 840.0087453600013, 843.2021915704545, 846.3956377809076, 849.5890839913608 $
              , 852.782530201814, 855.9759764122671, 859.1694226227203, 862.3628688331735, 865.5563150436266 $
              , 868.7497612540798, 871.943207464533, 875.1366536749862, 878.3300998854393, 881.5235460958925 $
              , 884.7169923063457, 887.9104385167989, 891.103884727252, 894.2973309377052, 897.4907771481583 $
              , 900.6842233586115, 903.8776695690647, 907.0711157795179, 910.264561989971, 913.4580082004242 $
              , 916.6514544108774, 919.8449006213305, 923.0383468317837, 926.2317930422369, 929.4252392526901 $
              , 932.6186854631432, 935.8121316735964, 939.0055778840496, 942.1990240945026, 945.3924703049559 $
              , 948.585916515409, 951.7793627258623, 952.8263378349274, 954.9728089363153, 958.1662551467685 $
              , 958.2751780026894, 961.3597013572218, 963.7240181704512, 964.5531475676748, 967.746593778128 $
              , 969.1728583382131, 970.9400399885812, 974.1334861990345, 974.6216985059751, 977.3269324094875 $
              , 980.0705386737371, 980.5203786199407, 983.713824830394, 985.5193788414991, 986.907271040847 $
              , 990.1007172513002, 990.968219009261, 993.2941634617534, 996.417059177023, 996.4876096722066 $
              , 1001.8658993447849, 1007.3147395125469, 1012.7635796803088, 1018.2124198480708, 1023.6612600158328 $
              , 1029.1101001835948, 1034.5589403513568, 1040.0077805191186, 1045.4566206868806, 1050.9054608546426 $
              , 1056.3543010224043, 1061.8031411901666, 1067.2519813579283, 1072.7008215256903, 1078.1496616934523 $
              , 1083.5985018612141, 1089.0473420289761, 1094.4961821967381, 1099.9450223645001, 1105.3938625322621 $
              , 1110.8427027000241, 1116.291542867786, 1121.740383035548, 1127.18922320331, 1132.638063371072 $
              , 1138.0869035388337, 1143.535743706596, 1148.9845838743577, 1154.4334240421197, 1159.8822642098817 $
              , 1165.3311043776434, 1170.7799445454057, 1176.2287847131674, 1181.6776248809294, 1187.1264650486914 $
              , 1192.5753052164534, 1198.0241453842152, 1203.4729855519772, 1208.9218257197392, 1214.3706658875012 $
              , 1219.819506055263, 1225.2683462230252, 1230.717186390787, 1236.166026558549, 1241.614866726311 $
              , 1247.0637068940728, 1252.512547061835, 1257.9613872295968, 1263.4102273973588, 1268.8590675651208 $
              , 1274.3079077328828, 1279.7567479006448, 1285.2055880684065, 1290.6544282361685, 1296.1032684039305 $
              , 1301.5521085716923, 1307.0009487394543, 1312.4497889072163, 1317.8986290749783, 1323.3474692427403 $
              , 1328.7963094105023, 1334.2451495782643, 1339.693989746026, 1345.142829913788, 1350.59167008155 $
              , 1356.040510249312, 1361.4893504170739, 1366.9381905848359, 1372.3870307525979, 1377.8358709203599 $
              , 1383.2847110881219, 1388.7335512558839, 1394.1823914236456, 1399.6312315914076, 1405.0800717591696 $
              , 1410.5289119269316, 1415.9777520946934, 1421.4265922624554, 1426.8754324302174, 1432.3242725979794 $
              , 1437.7731127657414, 1443.2219529335034, 1448.6707931012652, 1454.1196332690272, 1459.5684734367892 $
              , 1465.0173136045512, 1470.466153772313, 1475.914993940075, 1481.363834107837, 1486.812674275599 $
              , 1492.261514443361, 1497.7103546111227, 1503.1591947788847, 1508.6080349466467, 1514.0568751144087 $
              , 1519.5057152821707, 1524.9545554499325, 1530.4033956176945, 1535.8522357854565, 1541.3010759532185 $
              , 1546.7499161209805, 1552.1987562887423, 1557.6475964565043, 1563.0964366242663, 1568.5452767920283 $
              , 1573.9941169597903, 1579.442957127552, 1584.891797295314, 1590.340637463076, 1595.789477630838 $
              , 1601.2383177986, 1606.6871579663618, 1612.1359981341238, 1617.5848383018858, 1623.0336784696478 $
              , 1628.4825186374098, 1633.9313588051718, 1639.3801989729336, 1644.8290391406956, 1650.2778793084576 $
              , 1655.7267194762196, 1661.1755596439814, 1666.6243998117434, 1672.0732399795054, 1677.5220801472674 $
              , 1682.9709203150294, 1688.4197604827914, 1693.8686006505532, 1699.3174408183152, 1704.7662809860772 $
              , 1710.2151211538392, 1715.6639613216012, 1721.112801489363, 1726.561641657125, 1732.010481824887 $
              , 1737.459321992649, 1742.908162160411, 1748.3570023281727, 1753.8058424959347, 1759.2546826636967 $
              , 1764.7035228314587, 1770.1523629992207, 1775.6012031669825, 1781.0500433347445, 1786.4988835025065 $
              , 1791.9477236702685, 1797.3965638380305, 1802.8454040057923, 1808.2942441735543, 1813.7430843413163 $
              , 1819.1919245090783, 1824.6407646768403, 1830.089604844602, 1835.538445012364, 1840.987285180126 $
              , 1846.436125347888, 1851.88496551565, 1857.3338056834118, 1862.7826458511738, 1868.2314860189358 $
              , 1873.6803261866978, 1879.1291663544598, 1884.5780065222216, 1890.0268466899836, 1895.4756868577456 $
              , 1900.9245270255076, 1906.3733671932696, 1911.8222073610314, 1917.2710475287934, 1922.7198876965554 $
              , 1928.1687278643174, 1933.6175680320794, 1939.0664081998411, 1944.5152483676031, 1949.9640885353651 $
              , 1955.4129287031271, 1960.8617688708891, 1966.310609038651, 1971.759449206413, 1977.208289374175 $
              , 1982.657129541937, 1988.105969709699, 1993.5548098774607, 1999.0036500452227, 2004.4524902129847 $
              , 2009.9013303807467, 2015.3501705485087, 2020.7990107162705, 2026.2478508840325, 2031.6966910517945 $
              , 2037.1455312195565, 2042.5943713873182, 2048.0432115550802, 2053.492051722842, 2058.9408918906042 $
              , 2064.3897320583665, 2069.8385722261282, 2075.28741239389, 2080.7362525616522, 2086.185092729414 $
              , 2091.633932897176, 2097.082773064938, 2102.5316132327, 2107.9804534004616, 2113.429293568224 $
              , 2118.878133735986, 2124.326973903748, 2129.7758140715096, 2135.2246542392713, 2140.6734944070336 $
              , 2146.122334574796, 2151.5711747425576, 2157.0200149103193, 2162.4688550780816, 2167.9176952458433 $
              , 2173.366535413605, 2178.8153755813673, 2184.264215749129, 2189.713055916891, 2195.161896084653 $
              , 2200.6107362524153, 2206.059576420177, 2211.508416587939, 2216.9572567557007, 2222.406096923463 $
              , 2227.854937091225, 2233.303777258987, 2238.7526174267487, 2244.2014575945104, 2249.6502977622727 $
              , 2255.0991379300344, 2260.5479780977967, 2265.9968182655584, 2271.44565843332, 2276.8944986010824 $
              , 2282.343338768844, 2287.7921789366064, 2293.241019104368, 2298.68985927213, 2304.138699439892 $
              , 2309.587539607654, 2315.036379775416, 2320.485219943178, 2325.9340601109398, 2331.382900278702 $
              , 2336.8317404464638, 2342.280580614226, 2347.7294207819878, 2353.1782609497495, 2358.6271011175118 $
              , 2364.0759412852735, 2369.5247814530358, 2374.9736216207975, 2380.4224617885593, 2385.8713019563215 $
              , 2391.3201421240833, 2396.7689822918455, 2402.2178224596073, 2407.666662627369, 2413.1155027951313 $
              , 2418.564342962893, 2424.0131831306553, 2429.462023298417, 2434.910863466179, 2440.359703633941 $
              , 2445.808543801703, 2451.257383969465, 2456.706224137227, 2462.1550643049886, 2467.603904472751 $
              , 2473.0527446405126, 2478.501584808275, 2483.9504249760366, 2489.3992651437984, 2494.8481053115606 $
              , 2500.2969454793224, 2505.7457856470846, 2511.1946258148464, 2516.643465982608]
  return,wavelength
END

