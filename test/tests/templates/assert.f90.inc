#define PASTE(X)         X
#define PASTE2(X)        PASTE(X)_
#define CONCATHELP(X, Y) PASTE2(X)PASTE(Y)
#define CONCAT(X, Y)     CONCATHELP(X,Y)

#define ASSERT            CONCAT(assert,T)


!! assert x against y, z is true/false
function ASSERT(x,y) result(z)
  TT, intent(in) :: x,y
  logical             :: z

    z = (x .eq. y)

end function

#undef T
#undef TT

#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT

#undef ASSERT
