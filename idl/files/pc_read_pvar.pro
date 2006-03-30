; $Id: pc_read_pvar.pro,v 1.19 2006-03-30 13:50:27 ajohan Exp $
;
;   Read pvar.dat, or other PVAR file
;
pro pc_read_pvar, object=object, varfile=varfile, datadir=datadir, $
    npar_max=npar_max, quiet=quiet, qquiet=qquiet
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
;  Defaults.
;
default, datadir, 'data'
default, varfile, 'pvar.dat'
default, quiet, 0
default, qquiet, 0
if (qquiet) then quiet=1
;
;  Get necessary dimensions.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_read_pdim, obj=pdim, datadir=datadir, /quiet
pc_set_precision, dim=dim, /quiet
;
;  Derived dimensions.
;
mpvar=pdim.mpvar
npar=pdim.npar
default, npar_max, npar
ncpus=dim.nprocx*dim.nprocy*dim.nprocz
;
;  Time and grid parameters.
;
t=zero
x=fltarr(dim.mx)*one & y=fltarr(dim.my)*one & z=fltarr(dim.mz)*one
dx=zero &  dy=zero &  dz=zero
;
;
;
if (ncpus gt 1) then begin
  pc_read_dim, obj=procdim, datadir=datadir, proc=0, /quiet
endif else begin
  procdim=dim
endelse
xloc=fltarr(procdim.mx)*one & yloc=fltarr(procdim.my)*one & zloc=fltarr(procdim.mz)*one
;
;  Read variable indices from index.pro
;
default,datadir,'data'
cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),     "\n"'+"' "+datadir+'/index.pro'
spawn, cmd, result
res = flatten_strings(result)
if (execute(res) ne 1) then $
    message, 'There was a problem with index.pro', /INF
;
;  Define structure for data
;
varcontent=REPLICATE( $
    {varcontent_all_par, $
    variable   : 'UNKNOWN', $
    idlvar     : 'dummy', $
    idlinit    : 'fltarr(npar)*one', $
    skip       : 0}, $
    mpvar+1)

INIT_SCALAR  = 'fltarr(npar)*one'
INIT_3VECTOR = 'fltarr(npar,3)*one'
;
;  Go through all possible particle variables
;
varcontent[ixp].variable = 'Particle position (xx)'
varcontent[ixp].idlvar   = 'xx'
varcontent[ixp].idlinit  = INIT_3VECTOR
varcontent[ixp].skip     = 2

varcontent[ivpx].variable = 'Particle velocity (vv)'
varcontent[ivpx].idlvar   = 'vv'
varcontent[ivpx].idlinit  = INIT_3VECTOR
varcontent[ivpx].skip     = 2

varcontent[iap].variable = 'Particle radius (a)'
varcontent[iap].idlvar   = 'a'
varcontent[iap].idlinit  = INIT_SCALAR

varcontent[inptilde].variable = 'Particle internal number (nptilde)'
varcontent[inptilde].idlvar   = 'nptilde'
varcontent[inptilde].idlinit  = INIT_SCALAR

varcontent[0].variable    = 'UNKNOWN'
varcontent[0].idlvar      = 'UNKNOWN'
varcontent[0].idlinit     = '0.'
varcontent[0].skip        = 0
;
;  Put variable names in array
;
variables = (varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)[1:*]
;
;  Define arrays from contents of varcontent
;
totalvars = mpvar
for iv=1L,totalvars do begin
  if (varcontent[iv].variable eq 'UNKNOWN') then $
      message, 'Unknown variable at position ' + str(iv) $
      + ' needs declaring in pc_read_pvar.pro', /INF
  if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
      +' - '+ varcontent[iv].idlvar, /INFO
  iv=iv+varcontent[iv].skip
endfor
;
;  Define arrays for temporary storage of data.
;
array=fltarr(npar_max,totalvars)*one
ipar=lonarr(npar)
tarr=fltarr(ncpus)*one
t=zero
npar_loc=0L
;
;  Loop over processors.
;
for i=0,ncpus-1 do begin

  if (not keyword_set(quiet)) then $
      print,'Loading chunk ', strtrim(str(i+1)), ' of ', $
      strtrim(str(ncpus)), ' (', $
      strtrim(datadir+'/proc'+str(i)+'/'+varfile), ')...'

  filename=datadir+'/proc'+strtrim(i,2)+'/'+varfile 
;
;  Check if file exists.
;
  dummy=findfile(filename, COUNT=countfile)
  if (not countfile gt 0) then begin
    print, 'ERROR: cannot find file '+ filename
    stop
  endif
;
;  Get a unit number and open file.
;
  get_lun, file
  close, file
  openr, file, filename, /F77
;
;  Read the number of particles at the local processor together with their
;  global index numbers.
;
  readu, file, npar_loc 
;
;  Read particle data (if any).
;
  if (npar_loc ne 0) then begin
;    
    ipar_loc=lonarr(npar_loc)
    readu, file, ipar_loc
;
;  Register particle indices for later check if all particles have been read.
;  
    for k=0,npar_loc-1 do begin
      ipar[ipar_loc[k]-1]=ipar[ipar_loc[k]-1]+1
    endfor
;
;  Read local processor data.
;
    array_loc=fltarr(npar_loc,mpvar)*one
    readu, file, array_loc
;
;  Put local processor data into proper place in global data array
;        
    for k=0,npar_loc-1 do begin
      if (ipar_loc[k] lt npar_max) then array[ipar_loc[k]-1,*]=array_loc[k,*]
    endfor
;
  endif
;
;  Read time.
;
  readu, file, t, xloc, yloc, zloc, dx, dy, dz
;
;  Create global x, y and z arrays from local ones.
;
  if (ncpus gt 1) then begin
    pc_read_dim, object=procdim, datadir=datadir, proc=i, /quiet

    if (procdim.ipx eq 0L) then begin
      i0x=0L
      i1x=i0x+procdim.mx-1L
      i0xloc=0L
      i1xloc=procdim.mx-1L
    endif else begin
      i0x=procdim.ipx*procdim.nx+procdim.nghostx
      i1x=i0x+procdim.mx-1L-procdim.nghostx
      i0xloc=procdim.nghostx & i1xloc=procdim.mx-1L
    endelse
;
    if (procdim.ipy eq 0L) then begin
      i0y=0L
      i1y=i0y+procdim.my-1L
      i0yloc=0L
      i1yloc=procdim.my-1L
    endif else begin
      i0y=procdim.ipy*procdim.ny+procdim.nghosty
      i1y=i0y+procdim.my-1L-procdim.nghosty
      i0yloc=procdim.nghosty
      i1yloc=procdim.my-1L
    endelse
;
    if (procdim.ipz eq 0L) then begin
      i0z=0L
      i1z=i0z+procdim.mz-1L
      i0zloc=0L
      i1zloc=procdim.mz-1L
    endif else begin
      i0z=procdim.ipz*procdim.nz+procdim.nghostz
      i1z=i0z+procdim.mz-1L-procdim.nghostz
      i0zloc=procdim.nghostz
      i1zloc=procdim.mz-1L
    endelse

    x[i0x:i1x] = xloc[i0xloc:i1xloc]
    y[i0y:i1y] = yloc[i0yloc:i1yloc]
    z[i0z:i1z] = zloc[i0zloc:i1zloc]

  endif else begin

    x=xloc
    y=yloc
    z=zloc

  endelse
  tarr[i]=t
;
  close, file
  free_lun, file
;
endfor
;
;  Trim x, y and z arrays.
;
x=x[dim.l1:dim.l2]
y=y[dim.m1:dim.m2]
z=z[dim.n1:dim.n2]
;
;  Put data into sensibly named arrays.
;
for iv=1L,mpvar do begin
  res=varcontent[iv].idlvar+'=array[*,iv-1:iv-1+varcontent[iv].skip]'
  if (execute(res,0) ne 1) then $
    message, 'Error putting data into '+varcontent[iv].idlvar+' array'
  iv=iv+varcontent[iv].skip
endfor
;
;  Check if all particles found exactly once.
;
if ( (max(ipar) ne 1) or (min(ipar) ne 1)) then begin
  print, 'Warning: Some particles not found at all or found more'
  print, 'than once in snapshot files.'
  print, 'Particle number---No. of occurences'
  for i=0,npar-1 do begin
    if (ipar[i] ne 1) then begin
      print, i, ipar[i]
    endif
  endfor
endif
;
;  Put data and parameters in object.
;
makeobject="object = CREATE_STRUCT(name=objectname," + $
    "['t','x','y','z','dx','dy','dz'," + $
    arraytostring(variables,QUOTE="'",/noleader) + "]," + $
    "t,x,y,z,dx,dy,dz," + $
    arraytostring(variables,/noleader) + ")"
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine,object
endif

; If requested print a summary
;
;if keyword_set(STATS) or (not (keyword_set(NOSTATS) or keyword_set(quiet))) then begin
;  pc_object_stats,object,dim=dim,quiet=quiet
;endif
;
;  Check if times are consistent between processors.
;
for i=0,ncpus-1 do begin
  if (tarr[i] ne mean(tarr)) then begin
    print, 'The time of the snapshot at processor ', i, ' is strange!'
    print, 't[i], mean(t)=', tarr[i], mean(tarr)
  endif
endfor

if (not qquiet) then print,' t = ', mean(tarr)


end
