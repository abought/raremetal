cmake_minimum_required(VERSION 3.2)
project(dependencies)

include(ExternalProject)


ExternalProject_Add(libstatgen_project
                    BUILD_IN_SOURCE 1
                    GIT_REPOSITORY "https://github.com/statgen/libStatGen"
                    GIT_TAG "v1.0.14"
                    UPDATE_COMMAND ""
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND make
                    INSTALL_COMMAND "") #make "DESTDIR=${CMAKE_INSTALL_PREFIX}" install)

#ExternalProject_Add(libvc_project
#                    GIT_REPOSITORY "https://github.com/statgen/libvc"
#                    GIT_TAG "master"
#                    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")