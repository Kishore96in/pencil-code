;
;  $Id: del2.pro,v 1.2 2008-03-07 14:50:18 ajohan Exp $
;
;  Calculate div(grad(f)).
;
;  18-jan-08/anders: coded
;
function del2,f
  COMPILE_OPT IDL2,HIDDEN
;
  s=size(f)
;
  if (s[0] eq 3) then begin
;
    w=make_array(n_elements(f[*,0,0]),n_elements(f[0,*,0]),n_elements(f[0,0,*]),3,/nozero)
    w=xder2(f)+yder2(f)+zder2(f)
;
  endif else if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,/nozero)
    w=xder2(f)+yder2(f)+zder2(f)
;
  endif else begin
    print, 'error: del2 not implemented for arrays of size ', s
  endelse
;
  return, w
;
end
