#ifndef SIRE_VERSION_H
#define SIRE_VERSION_H

#define SIRE_REPOSITORY_URL       "@SVN_REPOSITORY_URL@"
#define SIRE_REPOSITORY_VERSION   "@SVN_VERSION_NUMBER@"
#define SIRE_REPOSITORY_BRANCH    "@SVN_REPOSITORY_BRANCH@"
#define SIRE_REPOSITORY_VERSION_IS_CLEAN "@SVN_IS_CLEAN@"

#define SIRE_VERSION_MAJOR        @S_VERSION_MAJOR@
#define SIRE_VERSION_MINOR        @S_VERSION_MINOR@
#define SIRE_VERSION_PATCH        @S_VERSION_PATCH@

#define SIRE_VERSION              @SIRE_VERSION_NUMBER@
#define SIRE_LIB_VERSION          @SIRE_VERSION_STRING@

#define SIRE_COMPILE_FLAGS        "@CMAKE_CXX_FLAGS@"
#define SIRE_LINK_FLAGS           "@CMAKE_SHARED_LINKER_FLAGS@"

#endif
