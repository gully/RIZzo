pro RIZzo_sandbox
;MGS

dir1='/Users/gully/Astronomy/BD/LIT_DATA/RIZzo/data/fits'
fn_list = 'fn_list.txt'
cd, dir1
dir4=dir1
readcol, fn_list, fns, format='A', delimiter=';'

n_files = n_elements(fns)

struct_template={$
fn:'x',$
TMASS_DESIG: 'x',$
SP_TYPE : 'x',$
SPT_NUM : 0.0,$
SPT_LCS : '_',$
SPT_REF : 'x',$
REF_URL : 'x'}
n=n_elements(fns)
s=replicate(struct_template, n)

help, s, /struct
s.fn=fns

for i=0, n-1 do begin
  d=readfits(s[i].fn, hdr, /silent)
  val=sxpar(hdr, 'SPT')
  print, val
  numspt=sptstr2num(val)
  s[i].SP_TYPE=val
  s[i].SPT_NUM=numspt
endfor

d=s

di= indgen(n)
ndw=n_elements(di)
;---------------------------------------------------------------------------
;***3) Make the list of spectral lines ***
H_ll=[6563.0, 8189.0, 7670.0, 7705.0] ;more lines!
r09_ll=H_ll
n_ll=n_elements(H_ll)

;Parameters for EW extraction determined by eye
dw_arr=[35.0,15.0, 10.0, 10.0]
delw_arr=[70.0,40.0,20.0, 20.0]
;---------------------------------------------------------------------------


;---------------------------------------------------------------------------
;***4) Set the parameters for the EW determination ***
cd, dir4
n_fs=n_elements(di)
ews=fltarr(n_fs, n_ll)
ewe=ews*0.0
;---------------------------------------------------------------------------

;---------------------------------------------------------------------------
;***5) Begin plotting functions and for loop ***  

for k=1, 2 do begin
;all the different objects
  for i=0, ndw-1 do begin
    dw=dw_arr[k]
    delw=delw_arr[k]
    
    ;dat=read_ascii(d[di[i]].fn, data_start=16, delimiter=' ')
    ;w=reform(dat.field1[0, *])
    ;f=reform(dat.field1[1, *])
    dat=readfits(d[i].fn)
    w=dat[*, 0]
    f=dat[*, 1]
    e=0.01*f
    n_p=n_elements(w)

    ;EW analysis
          ;get the nearby continuum indices, and then remove the avoid lines 
          ci_nb=where((w gt r09_ll[k]-dw-delw and w lt r09_ll[k]-dw) $
              or (w gt r09_ll[k]+dw and w lt r09_ll[k]+dw+delw) )
          ci=ci_nb;cmset_op(ci_nb, 'and', av_li, /not2)
          ri=where(w gt r09_ll[k]-dw-delw and w lt r09_ll[k]+dw+delw)
          ewi=where(w gt r09_ll[k]-dw and w lt r09_ll[k]+dw)
          if (ri[0] ne -1) and (ewi[0] ne -1) then begin
            
            ewiri=where(w[ri] gt r09_ll[k]-dw and w[ri] lt r09_ll[k]+dw)
            ak=linfit(w[ci], f[ci], covar=cv, yfit=cf_ci)
            Pk=[[1.0+fltarr(n_elements(ri))], [w[ri]]]
            cf=poly(w[ri], ak)
            nf=f[ri]/cf
            fci=f[ci]
            er=stddev(f[ci]/cf_ci)
            err=replicate(er, n_elements(ri))
            
            dli=(w[ewi+1]-w[ewi-1])/2.0
            ews[i, k]=total((1.0-nf[ewiri])*dli) ;from Cushing05
            ewe[i, k]=total( dli^2.0*(err[ewiri])^2.0 )
            ;ewe[i, k]=total( ((w[ewiri+1]-w[ewiri-1])/2.0)^2.0* ((err[ewiri]/f[ewiri])^2.0 )); + nf[ewiri]^2.0 * (sc/cf)^2.0 ) ) ;from Cushing05
          endif  
     endfor
endfor     
;---------------------------------------------------------------------------

;Write it out:
ew_out=ews[*, 2]+ews[*, 3]
ewe_out=(ewe[*, 2]^2.0+ewe[*, 3]^2.0)^0.5
forprint, d[di].fn, d[di].spt_num, ew_out, ewe_out, format='A,F,F,F',$
 textout='/Users/gully/Astronomy/BD/LIT_DATA/RIZzo/data/Rizzo_8189.dat'


d_spt_num=d.spt_num
d_ew=ew_out

device, decomposed=0
loadct, 13
!P.Multi = 0

outname='RIZzoPlot.eps'
psObject = Obj_New("FSC_PSConfig", /Color, /Times, /Bold, /Encapsulate, $
  Filename=outname, xsize=4.0, ysize=4.0)
thisDevice = !D.Name
Set_Plot, "PS"
Device, _Extra=psObject->GetKeywords()

;Colors
n1=10
c1=ceil(findgen(n1)/n1*312)
thick1=1.0


xtit='Spectral type (M0=0, L0=10)
ytit='EW (A)'
tit1='NaI 8189 A';'H'+greek('alpha', /append_font)+' (6563 A)'
 
plot, d_spt_num, d_ew, /nodata, title=tit1,  xtitle=xtit, ytitle=ytit,$
yrange=[-20, 20], xrange=[-5, 15], xstyle=1, ystyle=1
oplot, d_spt_num, d_ew, color=c1[9], psym=4, symsize=0.5  
       
;---------------------------------------------------------------------------
Device, /Close_File
Set_Plot, thisDevice
Obj_Destroy, psObject

print, 'END'


print, 'PAUSE'



end
