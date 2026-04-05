#pragma once

#include <volt/topology/crystal_coordination_topology.h>

#include <mutex>

namespace Volt{

inline void ensureCoordinationStructuresInitialized(){
    static std::once_flag initFlag;
    std::call_once(initFlag, []() {
        CoordinationStructures::initializeStructures();
    });
}

}
