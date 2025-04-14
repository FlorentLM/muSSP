cmake_minimum_required (VERSION 3.5)

project(muSSP)

set(SOURCES
    Graph.cpp
    main.cpp
    Node.cpp
    Sink.cpp
)

set(HEADERS
    Graph.h
    Node.h
    Sink.h
)

ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")