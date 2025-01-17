include(Util)
include(FetchContent)

FetchContent_Declare(
    LIBNEO
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG main
)
FetchContent_MakeAvailable(LIBNEO)

