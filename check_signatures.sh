#\!/bin/bash

echo "Checking function signatures..."

# Function enorm
echo -e "\n=== enorm ==="
echo "Interface:"
sed -n '/function enorm/,/end function/p' src/contrib/minpack_interfaces.f90  < /dev/null |  head -5
echo "Implementation:"
sed -n '/^function enorm/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine dogleg  
echo -e "\n=== dogleg ==="
echo "Interface:"
sed -n '/subroutine dogleg/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -6
echo "Implementation:"
sed -n '/^subroutine dogleg/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine fdjac1
echo -e "\n=== fdjac1 ==="
echo "Interface:"
sed -n '/subroutine fdjac1/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -7
echo "Implementation:"
sed -n '/^subroutine fdjac1/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine qrfac
echo -e "\n=== qrfac ==="
echo "Interface:"
sed -n '/subroutine qrfac/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -8
echo "Implementation:"
sed -n '/^subroutine qrfac/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine qform
echo -e "\n=== qform ==="
echo "Interface:"
sed -n '/subroutine qform/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -5
echo "Implementation:"
sed -n '/^subroutine qform/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine r1mpyq
echo -e "\n=== r1mpyq ==="
echo "Interface:"
sed -n '/subroutine r1mpyq/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -5
echo "Implementation:"
sed -n '/^subroutine r1mpyq/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine r1updt
echo -e "\n=== r1updt ==="
echo "Interface:"
sed -n '/subroutine r1updt/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -7
echo "Implementation:"
sed -n '/^subroutine r1updt/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

# Subroutine hybrd1
echo -e "\n=== hybrd1 ==="
echo "Interface:"
sed -n '/subroutine hybrd1/,/end subroutine/p' src/contrib/minpack_interfaces.f90 | head -8
echo "Implementation:"
sed -n '/^subroutine hybrd1/,/implicit none/{/implicit none/q;p}' src/contrib/minpack.f90

