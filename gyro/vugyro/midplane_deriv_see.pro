pro midplane_deriv_see, group=group

  common GLOBAL
  common MIDPLANE_DATA

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  if n_pwr eq 0 then return    
  if xregistered('midplane_deriv_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='plot', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='index up', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='index dn', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='+order', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='-order', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='lin/log', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=9)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
                  set_uvalue=state,$
                  /no_copy, $
                  /realize

  widget_control, draw, $
                  get_value=midplane_deriv_wid

  xmanager,'midplane_deriv_see', $
           base, $
           event='midplane_deriv_event',$
           group_leader=group

end
