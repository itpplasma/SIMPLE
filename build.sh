#!/ bin / bash -f


if [[ -d build ]]
then
rm - rf build
fi

if [[ $# -eq 0 ]]
then
export SIMPLE_COVERAGE = FALSE
else
    if [[ $1 == enable_tests ]]
    then
    export SIMPLE_COVERAGE = TRUE
    else
    echo "Type \"enable_test\" as first command line argument to enable pFUnit and code coverage."
    fi
fi


mkdir -p build
cd build
cmake ..
make -j



