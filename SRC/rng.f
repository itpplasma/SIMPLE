cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine random_num(ur)
      ur=zzg()
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine getran(irand,ur)
c
c  Produces the random number with zero deviation and unit square
c  deviation
c
c  Input parameters: irand - 0 for continious, 1 for +1 -1,
c  Output parameters: ur   - random number
c
      call random_num(ur)
c
      if(irand.eq.0) then
c
c  continiuos random number, constant is sqrt(12)
c
        ur=3.464102*(ur-.5)
      else
c
c  discrete random number
c
        if(ur.gt..5) then
          ur=1.
        else
          ur=-1.
        endif
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
