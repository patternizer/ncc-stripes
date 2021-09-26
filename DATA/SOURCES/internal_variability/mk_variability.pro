; For the Norfolk City Hall artwork stripes, we want a realisation of natural
; variability to superimpose on the simple climate model simulations of the
; future.

; The future projections come from the FAIR simple climate model, and this does
; not simulate internal variability -- and for the future projections it does
; not simulate variability driven by natural external forcings because these are
; not included in most future scenarios.  Therefore we want a plausible realisation
; of natural variability (internally generated and externally forced) to superimpose
; on the projected data.

; When we are dealing with global temperatures, we need natural variability for the 
; global mean.  So let's use the PAGES2k 2019 reconstruction of global mean temperature
; filtered to removed variability on timescales longer than 150 years.

; Read in PAGES2k 2019 data
; (copied from https://github.com/patternizer/ncc-stripes/blob/main/DATA/PAGES2k.txt)
; Read it as one string per line because they use 'NA' for missing, which makes it difficult
; to read in as numerica data.

fnin='pages2k_2019.txt'
print,'Reading data from: '+fnin

nyr=2000
nheader=5
nline=nheader+nyr
rawst=strarr(nline)
openr,lp,/get_lun,fnin
readf,lp,rawst
free_lun,lp

; extract data lines (skip header lines)
rawdat=rawst[nheader:*]

; process each line in turn for the 9 columns of data
; and store values into a data array, dealing with
; missing values appropriately

ncol=9
alldat=fltarr(ncol,nyr)*!values.f_nan

for iyr = 0 , nyr-1 do begin
 dat1=rawdat[iyr]
 datbits=strsplit(dat1,/extract)
 ml=where(datbits eq 'NA',nmiss)
 if nmiss gt 0 then datbits[ml]='-999.99'
 alldat[*,iyr]=float(datbits)
endfor

ml=where(alldat eq -999.99,nmiss)
if nmiss gt 0 then alldat[ml]=!values.f_nan

; extract year, instrumental data and median of the paleo ensemble

timey=reform(alldat[0,*])
instr=reform(alldat[1,*])
pag2k=reform(alldat[2,*])

; Prepare for plotting

ploton,'ncc_variability'
loadct,33
def_1color,!p.color<255,color='black'
def_1color,!p.background,color='white'
multi_plot,nrow=4,layout='large'
if !d.name eq 'X' then begin
  win,ysize=700
  !p.font=-1
endif else begin
  !p.font=0
  device,/helvetica,bold=0,font_size=14
endelse
makecolors,['magenta','red','deepblue','brown','orange','mdyellow','green',$
 'vlblue','mlblue','deepblue','lpurple','palepurple','black','mgrey'],indgen(14)+20

; plot the timeseries to see how they look

yr= [ min([pag2k,instr],/nan) , max([pag2k,instr],/nan) ]
plot,timey,pag2k,/nodata,$
 title='PAGES2k 2019 and Cowtan&Way global T',$
 xstyle=3,$
 ystyle=3,yrange=yr,ytitle='Temperature anomaly (!Uo!NC from 1961-1990)'

oplot,timey,pag2k,color=21
oplot,timey,instr,color=29

; high-pass filter the series with moderate truncation to reduce end effects,
; and plot them again

sm=150.
filter_cru2,sm,nan=0.3,tsin=instr,tshigh=instrh
filter_cru2,sm,nan=0.3,tsin=pag2k,tshigh=pag2kh

yr= [ min([pag2kh,instrh],/nan) , max([pag2kh,instrh],/nan) ]
plot,timey,pag2kh,/nodata,$
 title='PAGES2k 2019 and Cowtan&Way global T high-pass '+txt(sm)+' yr',$
 xstyle=3,$
 ystyle=3,yrange=yr,ytitle='Temperature anomaly (!Uo!NC)'

oplot,timey,pag2kh,color=21
oplot,timey,instrh,color=29

; inflate the variability of the Pages2k high-pass series so that its
; variance is closer to the instrumental high-pass variance

fac=1.6
pag2kh=pag2kh*fac

; calculate and plot SD of high-pass series in a running window

winlen=100
hwl=floor(winlen/2.)
nwin=nyr-winlen+1
instrhsd=fltarr(nyr)*!values.f_nan
pag2khsd=fltarr(nyr)*!values.f_nan
for iwin = 0 , nwin-1 do begin
 jwin=iwin+hwl    ; store in the mid-point of the window
 ts1=pag2kh[iwin:iwin+winlen-1]
 pag2khsd[jwin]=stddev(ts1)
 ts1=instrh[iwin:iwin+winlen-1]
 instrhsd[jwin]=stddev(ts1)
endfor

yr= [ min([pag2khsd,instrhsd],/nan) , max([pag2khsd,instrhsd],/nan) ]
plot,timey,pag2khsd,/nodata,$
 title='PAGES2k 2019 (inflated by '+txt(fac)+$
  ') and Cowtan&Way global T standard deviation of high-pass '+txt(sm)+' yr',$
 xstyle=3,$
 ystyle=3,yrange=yr,ytitle='Temperature SD (!Uo!NC)'

oplot,timey,pag2khsd,color=21
oplot,timey,instrhsd,color=29

; extract two 200-year segments of the inflated PAGES2k high-pass series, set the
; time-mean of each segment to zero, and output for use in superimposing on the
; FAIR simulations. also plot them

seglen=200
segst=[1350,1550]
nseg=n_elements(segst)

yr= [ min([pag2kh,instrh],/nan) , max([pag2kh,instrh],/nan) ]
plot,timey,pag2kh,/nodata,$
 title='PAGES2k 2019 (inflated by '+txt(fac)+') segments',$
 xstyle=3,$
 ystyle=3,yrange=yr,ytitle='Temperature anomaly (!Uo!NC)'

for iseg = 0 , nseg-1 do begin

 yr1=segst[iseg]
 iyr1=where(timey eq yr1) & iyr1=iyr1[0]

 x1=timey[iyr1:iyr1+seglen-1]
 ts1=pag2kh[iyr1:iyr1+seglen-1]
 mkanomaly,ts1
 oplot,x1,ts1,color=21+iseg

 datout=transpose([[x1],[ts1]])
 fnout='variability_realisation'+txt(iseg)+'.txt'
 print,'Writing variability realisation to: '+fnout
 openw,lv,/get_lun,fnout
 printf,lv,datout,format='(i6,f10.3)'
 free_lun,lv

endfor


; close plot and transfer to my website for viewing
;plotoff,/web
plotoff,/pdf,png=600,/notx

end
