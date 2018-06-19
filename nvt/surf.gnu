      u1(xx,yy) = 4.0*(xx**2+yy**2-1.0)**2*yy**2 #1.0
      u2(xx,yy) = -exp(-4.0*((xx-1)**2+yy**2)) # ..1.0
      u3(xx,yy) = -exp(-4.0*((xx+1)**2+yy**2)) #4.0 1.0
      u4(xx) = exp(8.0*(xx-1.5))
      u5(xx) = exp(-8.0*(xx+1.5))
      u6(yy) = exp(-4.0*(yy+.25)) #.25
      u7(xx) = 0.2*exp(-8.0*xx**2)
      pot(xx,yy) = u1(xx,yy)+u2(xx,yy)+u3(xx,yy)+u4(xx)+u5(xx)+u6(yy)+u7(xx)
set zr[-1:1.5]
set xr[-1.5:1.5]
set yr[-.8:1.6]
set isosamples 130,130
#set view map
#set palette defined (-1 '#000000', 10 '#ffffff')
splot pot(x,y) w pm3d
#,'lines.txt'u 1:2:3 w d lc 0
#replot'Bubbles/fort.88'u 1:2:(0.5) w p
#replot'Bubbles/fort.99'u 1:2:(0.5) w p ps 3 pt 7
#replot'Bubbles/fort.99'u 3:4:(0.5) w p ps 3 pt 7
#replot'Bubbles/path0'u 1:2:(0.5) w lp ps 3 pt 7
#replot'Bubbles/sup'u 1:2:(0.5) w lp ps 3 pt 7