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
                    INSTALL_COMMAND "")