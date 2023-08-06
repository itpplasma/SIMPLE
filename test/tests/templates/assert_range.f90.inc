#define PASTE(X)         X
#define PASTE2(X)        PASTE(X)_
#define CONCATHELP(X, Y) PASTE2(X)PASTE(Y)
#define CONCAT(X, Y)     CONCATHELP(X,Y)

#define ASSERT            CONCAT(assert,T)


!! assert x against y, z is true/false
function ASSERT(x,y,eps) result(z)
  TT, intent(in) :: x,y
  logical             :: z
  double precision :: eps
  
  if ((x < (y + eps)) .and. (x > (y - eps))) then
    z = .True.
  else
    z = .False.
  end if

end function

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef ASSERT
