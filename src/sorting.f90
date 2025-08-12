      subroutine sortin(a,ipoi,n)
!
      double precision a
!
      dimension a(n),ipoi(n)
!
      do 1 i=1,n
        ipoi(i)=i
 1    continue
!

      do 3 i=1,n-1
        do 2 j=i+1,n
          if(a(ipoi(j)).lt.a(ipoi(i))) then
            isave=ipoi(i)
            ipoi(i)=ipoi(j)
            ipoi(j)=isave
          endif
 2      continue 
 3    continue 
!   
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
